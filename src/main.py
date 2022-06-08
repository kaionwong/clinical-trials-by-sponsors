import sys
import os
import re
import time
from dateutil import parser
import datetime
from datetime import date, timedelta
import pandas as pd
import numpy as np
from pytrials.client import ClinicalTrials

base_dir = os.path.dirname(os.path.realpath(__file__))
processed_data_dir = base_dir + '\\data\\processed\\'
raw_data_dir = base_dir + '\\data\\raw\\'
output_dir = base_dir + '\\output\\'

# Example of fields - https://clinicaltrials.gov/api/gui/ref/crosswalks
'''fields=["NCTId", "Condition", "Keyword", "BriefTitle", "StatusVerifiedDate", "OverallStatus", "WhyStopped",
        "StartDate", "PrimaryCompletionDate", "CompletionDate", "LeadSponsorName", "IsFDARegulatedDrug",
        "IsFDARegulatedDevice", "IsUnapprovedDevice", "HasExpandedAccess", "ExpandedAccessNCTId",
        "DesignPrimaryPurpose", "Phase", "DesignAllocation", "EnrollmentCount", "EnrollmentType",
        "DesignObservationalModel", "PrimaryOutcomeMeasure", "EligibilityCriteria", "StudyPopulation",
        "LocationCountry", "LocationCity", "LocationStatus"]'''


def pandas_output_setting():
    """Set pandas output display setting"""
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', None)
    ##pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 170)
    pd.set_option('display.max_colwidth', None)
    pd.options.mode.chained_assignment = None  # default='warn'


def cleanup(txt_string):
    # 1) replace symbol with space; 2) remove leading and trailing space; 3) remove words from `delete_words`;
    # 4) replace single character with empty string"
    delete_words = ['infection', 'infections', 'usually', 'other', 'species', 'list of']
    txt_string = re.sub("[^a-zA-Z0-9-' \n\.]", ' ', txt_string)
    for sub in delete_words:
        txt_string = txt_string.replace(sub, ' ')
        txt_string = txt_string.replace('  ', ' ')

    txt_string = txt_string.lstrip()
    txt_string = txt_string.rstrip()
    txt_string = txt_string.replace('  ', ' ')
    txt_string = txt_string.replace('   ', ' ')
    txt_string = txt_string.replace('    ', ' ')

    if len(txt_string) >= 2:
        return txt_string


def create_infectious_disease_list() -> list:
    output_list = []
    remove_list = ['new', 'babies', 'children', 'mc', 'hand', 'foot', 'mono', 'multiple', None]
    infectious_disease_list_filename = 'infectious_disease_list.csv'
    input_dir = raw_data_dir + infectious_disease_list_filename
    df = pd.read_csv(input_dir)

    infectious_disease_agent_list = df['Infectious agent'].tolist()
    infectious_disease_common_name_list = df['Common name'].tolist()

    for string in infectious_disease_agent_list:
        string = string.lower()
        string_list = re.split('[,;()]', string)

        for i in string_list:
            i_sublist = re.split(' and ', i)

            for j in i_sublist:
                j = cleanup(j)
                output_list.append(j)

    for string in infectious_disease_common_name_list:
        string = string.lower()
        string_list = re.split('[,;()]', string)

        for i in string_list:
            i_sublist = re.split(' and ', i)

            for j in i_sublist:
                j = cleanup(j)
                output_list.append(j)

    # remove duplicates and unwanted entries
    output_list = list(dict.fromkeys(output_list))

    final_output_list = []
    for i in output_list:
        if i not in remove_list:
            final_output_list.append(i)

    return final_output_list


def create_condition_list_column(my_string):
    my_string_original = my_string.lower()
    exclusion_list = ['disease', 'diseases', 'prevention', 'healthy', 'volunteer', 'volunteers']
    output_list = []
    my_string = my_string.lower()

    string_list = re.split('[,;()]', my_string)

    for i in string_list:
        i_sublist = re.split(' ', i)

        for j in i_sublist:
            j = cleanup(j)
            if j and (j not in exclusion_list):
                output_list.append(j)

    output_list.append(my_string_original)
    return output_list


def compare_lists(col_list, comparison_list):
    if any(x in col_list for x in comparison_list):
        return True
    else:
        return False


def extract_year_from_date_string(date_string):
    return parser.parse(date_string).date().year


def create_pharma_company_list() -> list:
    pharma_company_list_filename_list = ['pharmaceutical_company_list.csv', 'large_pharma_company_list.csv']
    exclusion_list = ['dr.', 'Dr.', 'United', 'united', 'sun', 'Sun']
    exclude_from_first_name_entity_list = ['dr', 'Dr', 'dr.', 'Dr.', 'DR', 'DR.', 'AC', '180', 'BB', 'C4', 'Can', 'Day',
                                           'G1', 'IO', 'KV', 'Leap', 'Relay', 'X4', 'AIM', 'AN2', 'university',
                                           'University', 'Jiangsu', 'jiangsu']
    inclusion_list = ['united therapeutics', 'united biomedical', 'united medical specialties', 'jiangsu pacific']
    inclusion_list = [i.lower() for i in inclusion_list]
    output_list = []

    for pharma_company_list in pharma_company_list_filename_list:
        input_dir = raw_data_dir + pharma_company_list
        df = pd.read_csv(input_dir)

        company_list = df['Pharma Company Name'].tolist()

        for company_string in company_list:
            company_string_split = company_string.split(' ')
            if len(company_string_split) >= 2:
                if company_string_split[0] in exclude_from_first_name_entity_list:
                    output_list.append(company_string_split[1])  # add second name entity
                else:
                    output_list.append(company_string_split[0])  # add first name entity
            output_list.append(company_string)

    output_list = list(dict.fromkeys(output_list))  # remove duplicates
    final_list = [item for item in output_list if item not in exclusion_list]
    final_list = final_list + inclusion_list
    final_list.sort()

    return final_list


def extract_trial_info_from_api(api_obj, company_string, search_string, max_studies=1000,
                                infectious_disease_filter=False):
    """Warning: currently when the search results are more than 1000, it will only output 1000 (API limit).
    There may be a get-around by looping. However, for the time being, capping the starting year should help.
    And the interpretation shouldn't skew too much. Printed warning is in-place for alerting."""
    company_string = company_string.lower()
    study_count = api_obj.get_study_count(search_expr=search_string)
    if study_count > 1000:
        print('Warning: search string="{}" resulted in more than 1000 studies'.format(search_string))

    corona_fields = api_obj.get_study_fields(
        search_expr=search_string,
        fields=["NCTId", "LeadSponsorName", "Condition", "Keyword", "BriefTitle", "OverallStatus", "StatusVerifiedDate",
                "StartDate", "CompletionDate", "IsFDARegulatedDrug", "IsFDARegulatedDevice", "Phase", "LocationCountry",
                "PrimaryOutcomeMeasure"],
        max_studies=max_studies,
        fmt="csv",
    )

    df = pd.DataFrame.from_records(corona_fields[1:], columns=corona_fields[0])
    df['LeadSponsorNameLowercase'] = df['LeadSponsorName'].str.lower()
    df['CompanyAsLeadSponsor'] = df['LeadSponsorNameLowercase'].str.contains(company_string)

    # Match `Condition` with infectious disease list
    infectious_disease_list = create_infectious_disease_list()
    df['ConditionList'] = df['Condition'].apply(create_condition_list_column)
    df['LikelyIsInfectiousDisease'] = df['ConditionList'].apply(compare_lists, args=(infectious_disease_list,))

    # Final clean up
    df = df[df['CompanyAsLeadSponsor'] == True]
    df = df.drop(['LeadSponsorNameLowercase', 'CompanyAsLeadSponsor'], axis=1)

    if infectious_disease_filter:
        df = df[df['LikelyIsInfectiousDisease'] == True]

    total_count = len(df)
    df['TotalCount'] = total_count

    return df


def convert_string_date_into_date_format(string_date):
    """Expect string_date to be in the form: dd/mm/yyyy, such as 01/12/2020 or 28/02/2018"""
    return datetime.datetime.strptime(string_date, '%d/%m/%Y').date()


def convert_date_format_into_string_date(date):
    """Output date as string in the form: dd/mm/yyyy, such as 01/12/2020 or 28/02/2018"""
    return ('{}/{}/{}'.format(date.strftime('%d'), date.strftime('%m'), date.strftime('%Y')))


def main():
    # Search expression URL = https://clinicaltrials.gov/api/gui/ref/syntax
    pandas_output_setting()
    global_start_date = '01/01/2018'  # format dd/mm/yyyy (i.e., 28/02/2018)
    save_switch = True  # WARNING: will overwrite existing
    add_handcrafted_pharma_company_list_switch = True
    ct_obj = ClinicalTrials()
    max_studies_num = 1000  # 1000 is max per clinicaltrials.gov API call
    start_date = global_start_date
    day_segment = 60  # control how many days as 1 chunk in a search query; smaller chunk size, the longer it takes
    placeholder_company = '_placeholder_nonsensical_text'
    placeholder_search_string = '{} AND AREA[StartDate]RANGE[{}, MAX]'.format(placeholder_company, '01/01/2018')

    try:
        # If the file already exists
        df_master = pd.read_csv(output_dir + 'clinical_trials_by_pharma_companies_{}_onward.csv'.format(
            global_start_date.replace('/', '_')))

    except Exception as e:
        print(e)
        # Othersise, do first run to get the appropriate column labels onto an empty df
        df_master = extract_trial_info_from_api(api_obj=ct_obj, company_string=placeholder_company,
                                                search_string=placeholder_search_string,
                                                max_studies=max_studies_num,
                                                infectious_disease_filter=True)
        df_master = df_master[0:0]  # delete all rows but retain column names

    # Get pharma company list
    pharma_company_list = create_pharma_company_list()
    handcrafted_pharma_company_list = ['America Holdings', 'IQVIA', 'Syneos', 'Parexel', 'PRA',
                                       'Pharmaceutical Product Development', 'Charles River', 'Icon', 'WuXi Apptec',
                                       'Medpace', 'KCR', 'Pharm-Olam']
    handcrafted_pharma_company_list = [i.lower() for i in handcrafted_pharma_company_list]

    if add_handcrafted_pharma_company_list_switch:
        pharma_company_list += handcrafted_pharma_company_list

    # Loop through the list of pharma companies
    company_loop_count = 0
    for company in pharma_company_list[0:-1]:  # include all or a subset of companies
        try:
            # Create a fresh local df_combined which will be added to the df_master per company loop
            df_combined = extract_trial_info_from_api(api_obj=ct_obj, company_string=placeholder_company,
                                                      search_string=placeholder_search_string,
                                                      max_studies=max_studies_num,
                                                      infectious_disease_filter=True)
            df_combined = df_combined[0:0]  # delete all rows but retain column names

            df_excluded_companies = pd.read_csv(
                processed_data_dir + 'list_of_pharma_companies_without_infectious_disease_trial.csv')
            excluded_companies_list = df_excluded_companies['ExcludedCompanyName'].tolist()

            df_entered_companies = pd.read_csv(
                processed_data_dir + 'list_of_pharma_companies_entered_in_df_master.csv')
            entered_companies_list = df_entered_companies['EnteredCompanyName'].tolist()

            # If pharma companies already excluded (due to 0 clinical trials) or entered, skip to next loop item
            if (company in excluded_companies_list) | (company in entered_companies_list):
                continue

            company_loop_count += 1
            print('/Company loop count:{} /Current time:{} /Company:{}'.format(company_loop_count,
                                                                               datetime.datetime.now(),
                                                                               company))
            ct_obj = ClinicalTrials()
            full_start_date = convert_string_date_into_date_format(start_date)
            full_end_date = date.today() - timedelta(days=7)  # as a week before current day just to be safe
            num_day_difference = (full_end_date - full_start_date).days
            total_num_of_loops = int(np.ceil(num_day_difference / day_segment))
            local_start_date = full_start_date
            local_end_date = full_start_date + timedelta(days=day_segment)
            num_chunks_with_empty_row_count = 0

            for i in range(1, total_num_of_loops + 1):
                if i == total_num_of_loops:
                    curr_search_string = '{} AND AREA[StartDate]RANGE[{}, {}]'.format(company,
                                                                                      convert_date_format_into_string_date(
                                                                                          local_start_date),
                                                                                      convert_date_format_into_string_date(
                                                                                          full_end_date))
                    df = extract_trial_info_from_api(api_obj=ct_obj, company_string=company,
                                                     search_string=curr_search_string, max_studies=max_studies_num,
                                                     infectious_disease_filter=True)
                    df_combined = pd.concat([df_combined, df], ignore_index=True)

                    if len(df) == 0:
                        num_chunks_with_empty_row_count += 1

                    elif len(df) > 0:
                        if company not in entered_companies_list:
                            entered_companies_list.append(company)
                            entered_entry = pd.DataFrame([{'EnteredCompanyName': company}])
                            df_entered_companies = pd.concat([df_entered_companies, entered_entry], axis=0,
                                                             ignore_index=True)
                            df_entered_companies.to_csv(
                                processed_data_dir + 'list_of_pharma_companies_entered_in_df_master.csv',
                                index=False, header=True)

                    if num_chunks_with_empty_row_count == total_num_of_loops:
                        excluded_entry = pd.DataFrame([{'ExcludedCompanyName': company}])
                        df_excluded_companies = pd.concat([df_excluded_companies, excluded_entry], axis=0,
                                                          ignore_index=True)
                        df_excluded_companies.to_csv(
                            processed_data_dir + 'list_of_pharma_companies_without_infectious_disease_trial.csv',
                            index=False, header=True)

                    del df

                    break

                else:
                    curr_search_string = '{} AND AREA[StartDate]RANGE[{}, {}]'.format(company,
                                                                                      convert_date_format_into_string_date(
                                                                                          local_start_date),
                                                                                      convert_date_format_into_string_date(
                                                                                          local_end_date))
                    df = extract_trial_info_from_api(api_obj=ct_obj, company_string=company,
                                                     search_string=curr_search_string, max_studies=max_studies_num,
                                                     infectious_disease_filter=True)
                    df_combined = pd.concat([df_combined, df], ignore_index=True)

                    if len(df) == 0:
                        num_chunks_with_empty_row_count += 1

                    elif len(df) > 0:
                        if company not in entered_companies_list:
                            entered_companies_list.append(company)
                            entered_entry = pd.DataFrame([{'EnteredCompanyName': company}])
                            df_entered_companies = pd.concat([df_entered_companies, entered_entry], axis=0,
                                                             ignore_index=True)
                            df_entered_companies.to_csv(
                                processed_data_dir + 'list_of_pharma_companies_entered_in_df_master.csv',
                                index=False, header=True)

                    local_start_date += timedelta(days=day_segment)
                    local_end_date += timedelta(days=day_segment)

                    del df

                ##print(num_chunks_with_empty_row_count, total_num_of_loops)

            df_combined.drop_duplicates(subset=['LeadSponsorName', 'NCTId'], keep='first', ignore_index=True,
                                        inplace=True)
            print('/Row count of df_combined:{}'.format(
                len(df_combined)))  # check size to ensure if it doesn't exceed RAM
            df_combined['StartYear'] = df_combined['StartDate'].apply(extract_year_from_date_string)

        except Exception as e:
            wait_min = 30
            print(e, '/Waiting for {} min'.format(wait_min))
            time.sleep(wait_min * 60)

        if save_switch:
            try:
                df_master = pd.concat([df_master, df_combined], ignore_index=True)
                df_master.drop('Rank', inplace=True, axis=1)
                df_master.drop_duplicates(subset=['LeadSponsorName', 'NCTId'], keep='first', ignore_index=True,
                                          inplace=True)
                result_dir = output_dir + 'clinical_trials_by_pharma_companies_{}_onward.csv'.format(
                    global_start_date.replace('/', '_'))
                df_master.to_csv(result_dir, index=False, header=True)
            except Exception as e:
                print(e)


if __name__ == '__main__':
    main()

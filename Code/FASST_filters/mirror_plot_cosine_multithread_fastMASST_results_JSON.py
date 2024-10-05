import csv
import os
import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

def count_rows(filename):
    with open(filename, 'r') as file:
        return sum(1 for row in csv.reader(file, delimiter='\t')) - 1

def fetch_cosine_values(row, attempt=1):
    usi1 = row['USI_fastMASST']
    usi2 = row['USI_cluster']
    json_url = f"https://metabolomics-usi.gnps2.org/json/mirror/?usi1={usi1}&usi2={usi2}&width=10.0&height=6.0&mz_min=None&mz_max=None&max_intensity=125&annotate_precision=4&annotation_rotation=90&cosine=standard&fragment_mz_tolerance=0.1&grid=True&annotate_peaks=%5B%5B95.08549499511719%5D%2C%20%5B%5D%5D"
    interface_url = f"https://metabolomics-usi.gnps2.org/dashinterface/?usi1={usi1}&usi2={usi2}&width=10.0&height=6.0&mz_min=None&mz_max=None&max_intensity=125&annotate_precision=4&annotation_rotation=90&cosine=standard&fragment_mz_tolerance=0.1&grid=False"

    response = requests.get(json_url)
    if response.status_code == 200:
        json_data = response.json()
        cosine_value = json_data.get('cosine')
        precmz1 = json_data.get('spectrum1')['precursor_mz']
        precmz2 = json_data.get('spectrum2')['precursor_mz']
        return row, cosine_value, json_url, interface_url, True, precmz1, precmz2
    else:
        print(f"Attempt {attempt}: Error fetching data for standard cosine, usi1={usi1} and usi2={usi2}. Status code: {response.status_code}.", time.strftime('%H:%M:%S', time.localtime()))
        return row, None, json_url, interface_url, False

def load_processed_rows(output_tsv_file):
    processed = set()
    need_retry = []
    try:
        with open(output_tsv_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                if row['cosine_standard']:  # If there is a cosine value, add to processed
                    processed.add((row['USI_fastMASST'], row['USI_cluster']))
                else:  # If there is no cosine value, add to retry
                    need_retry.append(row)
    except FileNotFoundError:
        pass
    return processed, need_retry

def fetch_data_and_save_csv(input_tsv_file, output_tsv_file):
    processed_rows, rows_to_retry = load_processed_rows(output_tsv_file)
    print('Already processed:', len(processed_rows), 'Rows needing retry:', len(rows_to_retry))
    counter = 0

    # Check if file exists to determine if we need to write headers
    file_exists = os.path.exists(output_tsv_file)

    with open(input_tsv_file, 'r') as tsvfile, open(output_tsv_file, 'a', newline='') as outfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        fieldnames = ['group_annotation', 'USI_fastMASST', 'USI_cluster', 'url_standard', 'url_interface_standard', 'cosine_standard', 'precmz1', 'precmz2']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')

        # Only write the header if the file is new
        if not file_exists:
            writer.writeheader()

        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = {}
            # Resubmit rows needing retry
            for row in rows_to_retry:
                futures[executor.submit(fetch_cosine_values, row, attempt=2)] = row

            # Submit new rows for processing
            for row in reader:
                if (row['USI_fastMASST'], row['USI_cluster']) not in processed_rows:
                    futures[executor.submit(fetch_cosine_values, row)] = row

            for future in as_completed(futures):
                result = future.result()
                if result[4]:  # If fetch was successful
                    writer.writerow({
                        'group_annotation': result[0]['group_annotation'],
                        'USI_fastMASST': result[0]['USI_fastMASST'],
                        'USI_cluster': result[0]['USI_cluster'],
                        'url_standard': result[2],
                        'url_interface_standard': result[3],
                        'cosine_standard': result[1],
                        'precmz1': result[5],
                        'precmz2': result[6],
                    })
                    counter += 1
                    print(counter, "Current time:", time.strftime('%H:%M:%S', time.localtime()))
                    if counter % 20 == 0:
                        print(f"Processed {counter} requests so far.")

input_csv_file = '/.../Example_input.txt'
output_csv_file = '/.../output.tsv'
fetch_data_and_save_csv(input_csv_file, output_csv_file)
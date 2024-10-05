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
    svg_url = f'https://de.metabolomics-usi.gnps2.org/svg/mirror/?usi1={usi1}&usi2={usi2}'

    try:
        response = requests.get(svg_url, timeout=7)  # 5 second timeout
        if response.status_code == 200:
            cosine_value = response.text.split('Cosine similarity = ')[1].split('<')[0]
            try:
                precmz1 = response.text.split(r'Precursor $m$/$z$: ')[1].split(' ')[0]
            except:
                precmz1 = 0
            try:
                precmz2 = response.text.split(r'Precursor $m$/$z$: ')[2].split(' ')[0]
            except:
                precmz2 = 0
            return row, cosine_value, svg_url, True, precmz1, precmz2
        else:
            print(
                f"Attempt {attempt}: Error fetching data for SVG cosine, usi1={usi1} and usi2={usi2}. Status code: {response.status_code}")
            return row, None, svg_url, False
    except requests.exceptions.Timeout:
        print(f"Timeout occurred for usi1={usi1} and usi2={usi2}.")
        return row, None, svg_url, False
    except requests.exceptions.RequestException as e:
        print(f"Request exception for usi1={usi1} and usi2={usi2}: {str(e)}")
        return row, None, svg_url, False


def load_processed_rows(output_tsv_file):
    processed = set()
    need_retry = []
    try:
        with open(output_tsv_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                if row['cosine_SVG']:  # If there is a cosine value, add to processed
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

    # Check if file exists to decide whether to write headers
    file_exists = os.path.exists(output_tsv_file)

    with open(input_tsv_file, 'r') as tsvfile, open(output_tsv_file, 'a', newline='') as outfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        fieldnames = ['group_annotation', 'USI_fastMASST', 'USI_cluster', 'url_SVG', 'cosine_SVG', 'precmz1', 'precmz2']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')

        # Only write the header if the file is new
        if not file_exists:
            writer.writeheader()

        with ThreadPoolExecutor(max_workers=1) as executor:
            futures = {}
            for row in rows_to_retry:  # Retry previously failed rows
                futures[executor.submit(fetch_cosine_values, row)] = row

            for row in reader:  # Process new rows
                if (row['USI_fastMASST'], row['USI_cluster']) not in processed_rows:
                    futures[executor.submit(fetch_cosine_values, row)] = row

            for future in as_completed(futures):
                result = future.result()
                if result[3]:  # If fetch was successful
                    writer.writerow({
                        'group_annotation': result[0]['group_annotation'],
                        'USI_fastMASST': result[0]['USI_fastMASST'],
                        'USI_cluster': result[0]['USI_cluster'],
                        'url_SVG': result[2],
                        'cosine_SVG': result[1],
                        'precmz1': result[4],
                        'precmz2': result[5],
                    })
                    counter += 1
                    print(counter, time.strftime('%H:%M:%S', time.localtime()))
                    if counter % 1 == 0:
                        print(f"Processed {counter} requests so far.")

input_csv_file = '/.../Example_input.txt'
output_csv_file = '/.../output.tsv'
fetch_data_and_save_csv(input_csv_file, output_csv_file)

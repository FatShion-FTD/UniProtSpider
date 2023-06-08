# coding:utf-8
import csv
import requests
import logging
import json


logger = logging.getLogger()
logger.setLevel(logging.INFO) 
formatter = logging.Formatter('%(asctime)s :: %(filename)s :: %(funcName)s [line:%(lineno)d] :: %(levelname)s :: %(message)s')
sh = logging.StreamHandler()
sh.setFormatter(formatter)
logger.addHandler(sh)
fh = logging.FileHandler("script/uniprot-download_true_format_fasta_query__28CTCF_29-2023.06.07-17.36.58.14.log", encoding='utf8')
fh.setFormatter(formatter)
logger.addHandler(fh)


failed_ids = []


def get_function_content(id: str, key_words) -> int:
    """
    Get the function content of the protein from given id in JSON format
    """
    logging.info(f"Getting function content of {id}")
    url = "https://rest.uniprot.org/uniprotkb/" + id + ".json"
    response = requests.get(url)

    if response.status_code == 200:
        try:
            content = json.loads(response.text)
            function_content = content["comments"][0]["texts"][0]["value"].lower()
        except:
            logging.error(f"Failed to get function content of {id}")
            failed_ids.append(id)
            return -1

        # Check whether the function content contains the key words
        for keyword in key_words:
            if keyword in function_content:
                logging.info(f"Found keyword '{keyword}' in function content of {id}")
                return 1
            
        logging.info(f"No keyword found in function content of {id}")
        return 0
    else:
        logging.error(f"Failed to connect to {id}")
        failed_ids.append(id)
        return -1


def get_data_from_fasta(file_path):
    """
    Get the data from fasta file
    """
    # save id and sequence as list
    data = []
    cache = {}
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if not len(cache) == 0:
                    data.append(cache)
                cache = {}
                id = line.split("|")[1]         # split the line by '|'
                cache["id"] = id
                cache["description"] = line.replace("\n", "")
                cache["sequence"] = ""
            else:
                cache["sequence"] = cache["sequence"] + line.replace("\n", "")
    if not len(cache) == 0:
        data.append(cache)
    return data


if __name__ == "__main__":
    key_words = ["transcription factor", "transcription activator", "transcription repressor"]
    out_path = "script/CTCF.csv"
    in_path = "script/uniprot-download_true_format_fasta_query__28CTCF_29-2023.06.07-17.36.58.14.fasta"
    fail_path = "script/failed_ids.txt"

    entries = get_data_from_fasta(in_path)

    with open(out_path, "a", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["id", "description", "sequence", "TF/NTF"])
        writer.writeheader()
        for entry in entries:
            id = entry["id"]
            result = get_function_content(id, key_words)
            entry["TF/NTF"] = "TF" if result == 1 else "NTF"
            if result != -1:
                writer.writerow(entry)
            else:
                continue

    with open(fail_path, "w", newline='') as f:
        writer = csv.writer(f)
        for id in failed_ids:
            writer.writerow([id])
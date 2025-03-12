import pandas as pd
import ftplib
import os
import tarfile
import glob
import shutil
import time
import argparse

def construct_ftp_path(gse_id):
    num = int(gse_id[3:])
    nnn_group = f"GSE{num//1000}nnn"
    return f"/geo/series/{nnn_group}/{gse_id}/suppl/{gse_id}_RAW.tar"

def extract_and_filter_idats(gse_id, samples, local_dir):
    tar_file = os.path.join(local_dir, f"{gse_id}_RAW.tar")
    temp_dir = os.path.join(local_dir, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    
    # Extract tar file
    with tarfile.open(tar_file) as tar:
        tar.extractall(temp_dir)
    
    # Find all .idat or .idat.gz files
    idat_files = []
    for ext in ['*.idat', '*.idat.gz']:
        idat_files.extend(glob.glob(os.path.join(temp_dir, '**', ext), recursive=True))
    
    # Filter for specific samples
    kept_files = []
    for sample in samples:
        matching_files = [f for f in idat_files if sample in f]
        kept_files.extend(matching_files)
    
    # Move matching files to final location
    for file in kept_files:
        shutil.move(file, os.path.join(local_dir, os.path.basename(file)))
    
    # Cleanup
    shutil.rmtree(temp_dir)
    os.remove(tar_file)

def download_gse_datasets():
    parser = argparse.ArgumentParser()
    parser.add_argument('--annotations', type=str, required=True, help='Path to annotations CSV file')
    args = parser.parse_args()
    
    df = pd.read_csv(args.annotations)
    gse_groups = df.groupby('Dataset')['Sample'].apply(list).to_dict()
    
    ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    
    # Keep track of last successful GSE
    last_successful_gse = None
    
    for gse_id, samples in gse_groups.items():
        if not gse_id.startswith('GSE'):
            continue
        
        local_dir = f"../raw/GEO_zeroshot/{gse_id}"
        if os.path.exists(local_dir):
            print(f"Skipping {gse_id} as it already exists.")
            last_successful_gse = gse_id
            continue

        # Retry mechanism with max 3 retries
        retries = 3
        success = False
        while not success and retries > 0:
            try:
                ftp_path = construct_ftp_path(gse_id)
                os.makedirs(local_dir, exist_ok=True)
                local_file = os.path.join(local_dir, f"{gse_id}_RAW.tar")
                
                print(f"Downloading {gse_id}...")
                with open(local_file, 'wb') as fp:
                    ftp.retrbinary(f'RETR {ftp_path}', fp.write)
                
                print(f"Extracting and filtering files for {gse_id}...")
                extract_and_filter_idats(gse_id, samples, local_dir)
                
                success = True
                last_successful_gse = gse_id
            except Exception as e:
                retries -= 1
                if retries > 0:
                    print(f"Error downloading {gse_id}: {e}. Retrying in 10 seconds... ({retries} retries left)")
                    time.sleep(10)  # Retry after 10 seconds
                else:
                    print(f"Error downloading {gse_id}: {e}. Max retries reached. Skipping.")
                    if os.path.exists(local_file):
                        os.remove(local_file)
                    if os.path.isdir(local_dir) and not os.listdir(local_dir):
                        os.rmdir(local_dir)
                    success = True  # Skip this GSE after max retries

    print(f"Completed processing. Last successful GSE: {last_successful_gse}")

if __name__ == "__main__":
    download_gse_datasets()

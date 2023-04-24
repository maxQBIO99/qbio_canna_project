#!/bin/bash
# This script downloads the SRR files from a self-made textfile containing all the necessary links.
# After downloading the files it will compare the corresponding md5sum which are contained in another self-made textfile.

urls_file="urls.txt"
md5sums_file="md5sums_ena.txt"

# Read the URLs and md5sums into arrays
mapfile -t urls < "$urls_file"
mapfile -t md5sums < "$md5sums_file"

total_files=${#urls[@]}

for i in "${!urls[@]}"; do
  url="${urls[$i]}"
  expected_md5="${md5sums[$i]}"
  filename=$(basename "$url")

  echo "Processing file $(($i + 1)) of $total_files: $filename"

  # Download the file if it doesn't exist or if the md5sum is incorrect
  while true; do
    if [ ! -f "$filename" ]; then
      echo "Downloading $filename..."
      wget "$url" -O "$filename"
    fi

    # Check the md5sum of the downloaded file
    actual_md5=$(md5sum "$filename" | awk '{print $1}')

    if [ "$actual_md5" != "$expected_md5" ]; then
      echo "md5sum mismatch for $filename, redownloading..."
      rm "$filename"
    else
      echo "md5sum verified for $filename"
      break
    fi
  done

  echo "Done with file $(($i + 1)) of $total_files"
  echo
done

echo "All files downloaded and verified."

#!/bin/bash

set -e
for dir in Logs PDFs Images; do
  if [ -d "$dir" ]; then
    echo "Directory $dir already exists, skipping creation."
  else
     mkdir "$dir"
     chmod -R ug=rwx,o=rx "$dir"
     echo "Directory $dir created successfully."
  fi
done

run_workbook_template() {
  local template_name=$1
  local template_location=$2
  local template_count=$3
  echo "############ ${template_location}/${template_count} ${template_name}.R Starting... ########"
  Rscript ${template_name}.R 2>&1 | tee ./Logs/${template_name}.log; exit_status=${PIPESTATUS[0]}
  if [ $exit_status -ne 0 ]; then
    exit $exit_status
  fi
  if [ -f Rplots.pdf ]; then
    mv Rplots.pdf ./PDFs/${template_name}.pdf
  fi
  find . -maxdepth 1 -type f \( -iname '*.png' -o -iname '*.jpg' \) -exec mv {} Images/ \;
  echo "############ ${template_location}/${template_count} ${template_name}.R Completed! ########"
}
############ 1/24 template_PostIt.R ########
run_workbook_template "template_PostIt" "1" "24"

############ 2/24 template_PostIt_CellTypes.R ########
run_workbook_template "template_PostIt_CellTypes" "2" "24"

############ 3/24 template_PostIt_CombNorm.R ########
run_workbook_template "template_PostIt_CombNorm" "3" "24"

############ 4/24 template_PostIt_DualLabel.R ########
run_workbook_template "template_PostIt_DualLabel" "4" "24"

############ 5/24 template_PostIt_ProcessedInput.R ########
run_workbook_template "template_PostIt_ProcessedInput" "5" "24"

############ 6/24 template_PostIt_QCFiltered.R ########
run_workbook_template "template_PostIt_QCFiltered" "6" "24"

############ 7/24 template_PostIt_YourDataHere.R ########
run_workbook_template "template_PostIt_YourDataHere" "7" "24"

############ 8/24 template_ProcessedInputSO.R ########
run_workbook_template "template_ProcessedInputSO" "8" "24"

############ 9/24 template_QCFilteredSO.R ########
run_workbook_template "template_QCFilteredSO" "9" "24"

############ 10/24 template_CombNormSO.R ########
run_workbook_template "template_CombNormSO" "10" "24"

############ 11/24 template_MetadataTable_CombNorm.R ########
run_workbook_template "template_MetadataTable_CombNorm" "11" "24"

############ 12/24 template_CellTypesSO.R ########
run_workbook_template "template_CellTypesSO" "12" "24"

############ 13/24 template_MetadataTable_CellTypes.R ########
run_workbook_template "template_MetadataTable_CellTypes" "13" "24"

############ 14/24 template_SampleNames_CellTypes.R ########
run_workbook_template "template_SampleNames_CellTypes" "14" "24"

############ 15/24 template_DualLabelSO_CD79bposCD3neg_Bcells.R ########
run_workbook_template "template_DualLabelSO_CD79bposCD3neg_Bcells" "15" "24"

############ 16/24 template_FilteredSO.R ########
run_workbook_template "template_FilteredSO" "16" "24"

############ 17/24 template_MDT_DL_CD79bposCD3neg_Bcells.R ########
run_workbook_template "template_MDT_DL_CD79bposCD3neg_Bcells" "17" "24"

############ 18/24 template_FSO_CD79bposCD3neg_Bcells.R ########
run_workbook_template "template_FSO_CD79bposCD3neg_Bcells" "18" "24"

############ 19/24 template_MDT_FSO_CD79bposCD3neg_Bcells.R ########
run_workbook_template "template_MDT_FSO_CD79bposCD3neg_Bcells" "19" "24"

############ 20/24 template_Reclustered_CD79bposCD3neg_Bcells.R ########
run_workbook_template "template_Reclustered_CD79bposCD3neg_Bcells" "20" "24"

############ 21/24 template_RMDT_CD79bposCD3neg_Bcells.R ########
run_workbook_template "template_RMDT_CD79bposCD3neg_Bcells" "21" "24"

############ 22/24 template_RSN_CD79bposCD3neg_Bcells.R ########
run_workbook_template "template_RSN_CD79bposCD3neg_Bcells" "22" "24"

############ 23/24 template_DEG_CD79bposCD3neg_Bcells.R ########
run_workbook_template "template_DEG_CD79bposCD3neg_Bcells" "23" "24"

############ 24/24 template_SuppFIGURE_5E.R ########
run_workbook_template "template_SuppFIGURE_5E" "24" "24"


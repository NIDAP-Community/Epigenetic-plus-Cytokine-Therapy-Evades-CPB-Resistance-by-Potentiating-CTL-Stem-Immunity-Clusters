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
############ 1/42 template_PostIt_CellTypes.R ########
run_workbook_template "template_PostIt_CellTypes" "1" "42"

############ 2/42 template_PostIt_ColorByGene.R ########
run_workbook_template "template_PostIt_ColorByGene" "2" "42"

############ 3/42 template_PostIt_ColorByGeneList.R ########
run_workbook_template "template_PostIt_ColorByGeneList" "3" "42"

############ 4/42 template_PostIt_CombNorm.R ########
run_workbook_template "template_PostIt_CombNorm" "4" "42"

############ 5/42 template_PostIt_Filtered.R ########
run_workbook_template "template_PostIt_Filtered" "5" "42"

############ 6/42 template_PostIt_GenesClusters.R ########
run_workbook_template "template_PostIt_GenesClusters" "6" "42"

############ 7/42 template_PostIt_ModScore.R ########
run_workbook_template "template_PostIt_ModScore" "7" "42"

############ 8/42 template_PostIt_NameClusters.R ########
run_workbook_template "template_PostIt_NameClusters" "8" "42"

############ 9/42 template_PostIt_ProcessedInput.R ########
run_workbook_template "template_PostIt_ProcessedInput" "9" "42"

############ 10/42 template_PostIt_QCFiltered.R ########
run_workbook_template "template_PostIt_QCFiltered" "10" "42"

############ 11/42 template_PostIt_Reclustered.R ########
run_workbook_template "template_PostIt_Reclustered" "11" "42"

############ 12/42 template_PostIt_YourDataHere.R ########
run_workbook_template "template_PostIt_YourDataHere" "12" "42"

############ 13/42 template_ProcessedInputSO.R ########
run_workbook_template "template_ProcessedInputSO" "13" "42"

############ 14/42 template_QCFilteredSO.R ########
run_workbook_template "template_QCFilteredSO" "14" "42"

############ 15/42 template_CombNormSO.R ########
run_workbook_template "template_CombNormSO" "15" "42"

############ 16/42 template_MetadataTable_CombNorm.R ########
run_workbook_template "template_MetadataTable_CombNorm" "16" "42"

############ 17/42 template_CellTypesSO.R ########
run_workbook_template "template_CellTypesSO" "17" "42"

############ 18/42 template_MetadataTable_CellTypes.R ########
run_workbook_template "template_MetadataTable_CellTypes" "18" "42"

############ 19/42 template_ModScoreSO.R ########
run_workbook_template "template_ModScoreSO" "19" "42"

############ 20/42 template_SampleNames_CellTypes.R ########
run_workbook_template "template_SampleNames_CellTypes" "20" "42"

############ 21/42 template_SuppFIGURE_5Bpart1.R ########
run_workbook_template "template_SuppFIGURE_5Bpart1" "21" "42"

############ 22/42 template_SuppFIGURE_5Bpart2.R ########
run_workbook_template "template_SuppFIGURE_5Bpart2" "22" "42"

############ 23/42 template_SuppFIGURE_5Bpart3.R ########
run_workbook_template "template_SuppFIGURE_5Bpart3" "23" "42"

############ 24/42 template_ColorByGeneList_Old.R ########
run_workbook_template "template_ColorByGeneList_Old" "24" "42"

############ 25/42 template_FilteredSO.R ########
run_workbook_template "template_FilteredSO" "25" "42"

############ 26/42 template_MetadataTable_Filtered_9_10.R ########
run_workbook_template "template_MetadataTable_Filtered_9_10" "26" "42"

############ 27/42 template_MetadataTable_ModScore.R ########
run_workbook_template "template_MetadataTable_ModScore" "27" "42"

############ 28/42 template_NameClustersSO.R ########
run_workbook_template "template_NameClustersSO" "28" "42"

############ 29/42 template_ReclusteredSO.R ########
run_workbook_template "template_ReclusteredSO" "29" "42"

############ 30/42 template_SampleNames_Reclustered.R ########
run_workbook_template "template_SampleNames_Reclustered" "30" "42"

############ 31/42 template_SuppFIGURE_2E_ALL.R ########
run_workbook_template "template_SuppFIGURE_2E_ALL" "31" "42"

############ 32/42 template_SuppFIGURE_2E_ENT.R ########
run_workbook_template "template_SuppFIGURE_2E_ENT" "32" "42"

############ 33/42 template_SuppFIGURE_2E_ENTN803.R ########
run_workbook_template "template_SuppFIGURE_2E_ENTN803" "33" "42"

############ 34/42 template_SuppFIGURE_2E_ENTN803aPD1.R ########
run_workbook_template "template_SuppFIGURE_2E_ENTN803aPD1" "34" "42"

############ 35/42 template_SuppFIGURE_2E_ENTaPD1.R ########
run_workbook_template "template_SuppFIGURE_2E_ENTaPD1" "35" "42"

############ 36/42 template_SuppFIGURE_2E_N803.R ########
run_workbook_template "template_SuppFIGURE_2E_N803" "36" "42"

############ 37/42 template_SuppFIGURE_2E_N803aPD1.R ########
run_workbook_template "template_SuppFIGURE_2E_N803aPD1" "37" "42"

############ 38/42 template_SuppFIGURE_2E_PBS.R ########
run_workbook_template "template_SuppFIGURE_2E_PBS" "38" "42"

############ 39/42 template_SuppFIGURE_2E_aPD1.R ########
run_workbook_template "template_SuppFIGURE_2E_aPD1" "39" "42"

############ 40/42 template_MetadataTable_NameClusters.R ########
run_workbook_template "template_MetadataTable_NameClusters" "40" "42"

############ 41/42 template_MetadataTable_Reclustered.R ########
run_workbook_template "template_MetadataTable_Reclustered" "41" "42"

############ 42/42 template_DotGenesClusters.R ########
run_workbook_template "template_DotGenesClusters" "42" "42"


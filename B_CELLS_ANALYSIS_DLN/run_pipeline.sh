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
############ 1/31 template_PostIt_CellTypes.R ########
run_workbook_template "template_PostIt_CellTypes" "1" "31"

############ 2/31 template_PostIt_CombNorm.R ########
run_workbook_template "template_PostIt_CombNorm" "2" "31"

############ 3/31 template_PostIt_DualLabel.R ########
run_workbook_template "template_PostIt_DualLabel" "3" "31"

############ 4/31 template_PostIt_ProcessedInput.R ########
run_workbook_template "template_PostIt_ProcessedInput" "4" "31"

############ 5/31 template_PostIt_QCFiltered.R ########
run_workbook_template "template_PostIt_QCFiltered" "5" "31"

############ 6/31 template_PostIt_YourDataHere.R ########
run_workbook_template "template_PostIt_YourDataHere" "6" "31"

############ 7/31 template_ProcessedInputSO.R ########
run_workbook_template "template_ProcessedInputSO" "7" "31"

############ 8/31 template_QCFilteredSO.R ########
run_workbook_template "template_QCFilteredSO" "8" "31"

############ 9/31 template_CombNormSO.R ########
run_workbook_template "template_CombNormSO" "9" "31"

############ 10/31 template_MetadataTable_CombNorm.R ########
run_workbook_template "template_MetadataTable_CombNorm" "10" "31"

############ 11/31 template_CellTypesSO.R ########
run_workbook_template "template_CellTypesSO" "11" "31"

############ 12/31 template_MetadataTable_CellTypes.R ########
run_workbook_template "template_MetadataTable_CellTypes" "12" "31"

############ 13/31 template_SampleNames_CellTypes.R ########
run_workbook_template "template_SampleNames_CellTypes" "13" "31"

############ 14/31 template_DP_Bcells.R ########
run_workbook_template "template_DP_Bcells" "14" "31"

############ 15/31 template_DP_Bcells_MDT.R ########
run_workbook_template "template_DP_Bcells_MDT" "15" "31"

############ 16/31 template_DualLabelSO_Cd19posCd3neg.R ########
run_workbook_template "template_DualLabelSO_Cd19posCd3neg" "16" "31"

############ 17/31 template_Filteredso_B_Cells.R ########
run_workbook_template "template_Filteredso_B_Cells" "17" "31"

############ 18/31 template_MetadataTable_Cd19posCd3neg.R ########
run_workbook_template "template_MetadataTable_Cd19posCd3neg" "18" "31"

############ 19/31 template_DP_Bcells_FilteredSO.R ########
run_workbook_template "template_DP_Bcells_FilteredSO" "19" "31"

############ 20/31 template_FilteredDPB_MDT.R ########
run_workbook_template "template_FilteredDPB_MDT" "20" "31"

############ 21/31 template_FilteredSO_Cd19posCd3neg.R ########
run_workbook_template "template_FilteredSO_Cd19posCd3neg" "21" "31"

############ 22/31 template_MetadataTablefiltered_Cd19posCd3neg.R ########
run_workbook_template "template_MetadataTablefiltered_Cd19posCd3neg" "22" "31"

############ 23/31 template_ReclusterSO1.R ########
run_workbook_template "template_ReclusterSO1" "23" "31"

############ 24/31 template_ReclusterSO_Cd19posCd3neg.R ########
run_workbook_template "template_ReclusterSO_Cd19posCd3neg" "24" "31"

############ 25/31 template_SampleNames_Cd19posCd3neg.R ########
run_workbook_template "template_SampleNames_Cd19posCd3neg" "25" "31"

############ 26/31 template_RMDT_Cd19posCd3neg.R ########
run_workbook_template "template_RMDT_Cd19posCd3neg" "26" "31"

############ 27/31 template_RSO_MDT_BcellsDP.R ########
run_workbook_template "template_RSO_MDT_BcellsDP" "27" "31"

############ 28/31 template_RSO_SN_BcellsDP.R ########
run_workbook_template "template_RSO_SN_BcellsDP" "28" "31"

############ 29/31 template_DEGMarkers_Cd19posCd3neg.R ########
run_workbook_template "template_DEGMarkers_Cd19posCd3neg" "29" "31"

############ 30/31 template_DEG_BcellsDP.R ########
run_workbook_template "template_DEG_BcellsDP" "30" "31"

############ 31/31 template_SuppFIGURE_5D.R ########
run_workbook_template "template_SuppFIGURE_5D" "31" "31"


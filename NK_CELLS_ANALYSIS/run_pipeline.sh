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
############ 1/42 template_PostIt.R ########
run_workbook_template "template_PostIt" "1" "42"

############ 2/42 template_PostIt_1.R ########
run_workbook_template "template_PostIt_1" "2" "42"

############ 3/42 template_PostIt_CellTypes.R ########
run_workbook_template "template_PostIt_CellTypes" "3" "42"

############ 4/42 template_PostIt_CombNorm.R ########
run_workbook_template "template_PostIt_CombNorm" "4" "42"

############ 5/42 template_PostIt_DualLabel.R ########
run_workbook_template "template_PostIt_DualLabel" "5" "42"

############ 6/42 template_PostIt_ProcessedInput.R ########
run_workbook_template "template_PostIt_ProcessedInput" "6" "42"

############ 7/42 template_PostIt_QCFiltered.R ########
run_workbook_template "template_PostIt_QCFiltered" "7" "42"

############ 8/42 template_PostIt_YourDataHere.R ########
run_workbook_template "template_PostIt_YourDataHere" "8" "42"

############ 9/42 template_ProcessedInputSO.R ########
run_workbook_template "template_ProcessedInputSO" "9" "42"

############ 10/42 template_QCFilteredSO.R ########
run_workbook_template "template_QCFilteredSO" "10" "42"

############ 11/42 template_CombNormSO.R ########
run_workbook_template "template_CombNormSO" "11" "42"

############ 12/42 template_MetadataTable_CombNorm.R ########
run_workbook_template "template_MetadataTable_CombNorm" "12" "42"

############ 13/42 template_CellTypesSO.R ########
run_workbook_template "template_CellTypesSO" "13" "42"

############ 14/42 template_MetadataTable_CellTypes.R ########
run_workbook_template "template_MetadataTable_CellTypes" "14" "42"

############ 15/42 template_SampleNames_CellTypes.R ########
run_workbook_template "template_SampleNames_CellTypes" "15" "42"

############ 16/42 template_ColorByMeta.R ########
run_workbook_template "template_ColorByMeta" "16" "42"

############ 17/42 template_DualLabelSO_CD49bNkp46DP.R ########
run_workbook_template "template_DualLabelSO_CD49bNkp46DP" "17" "42"

############ 18/42 template_DLMDT_CD49bNkp46DP.R ########
run_workbook_template "template_DLMDT_CD49bNkp46DP" "18" "42"

############ 19/42 template_FSO_CD49bNkp46DP.R ########
run_workbook_template "template_FSO_CD49bNkp46DP" "19" "42"

############ 20/42 template_FSOMDT_CD49bNkp46DP.R ########
run_workbook_template "template_FSOMDT_CD49bNkp46DP" "20" "42"

############ 21/42 template_RSO_CD49bNkp46DP.R ########
run_workbook_template "template_RSO_CD49bNkp46DP" "21" "42"

############ 22/42 template_RMDT_CD49bNkp46DP.R ########
run_workbook_template "template_RMDT_CD49bNkp46DP" "22" "42"

############ 23/42 template_RSN_CD49bNkp46DP.R ########
run_workbook_template "template_RSN_CD49bNkp46DP" "23" "42"

############ 24/42 template_DEGMarkers_CD49bNkp46DP.R ########
run_workbook_template "template_DEGMarkers_CD49bNkp46DP" "24" "42"

############ 25/42 template_FIGURE_4F.R ########
run_workbook_template "template_FIGURE_4F" "25" "42"

############ 26/42 template_Il1bnegS100a9neg.R ########
run_workbook_template "template_Il1bnegS100a9neg" "26" "42"

############ 27/42 template_Metadata_Il1bnegS100a9neg.R ########
run_workbook_template "template_Metadata_Il1bnegS100a9neg" "27" "42"

############ 28/42 template_SuppFIGURE_5A_TripletvsENTN803.R ########
run_workbook_template "template_SuppFIGURE_5A_TripletvsENTN803" "28" "42"

############ 29/42 template_SuppFIGURE_5A_TripletvsN803.R ########
run_workbook_template "template_SuppFIGURE_5A_TripletvsN803" "29" "42"

############ 30/42 template_SuppFIGURE_5A_TripletvsN803aPD1.R ########
run_workbook_template "template_SuppFIGURE_5A_TripletvsN803aPD1" "30" "42"

############ 31/42 template_FilteredSO_Il1bnegS100a9neg.R ########
run_workbook_template "template_FilteredSO_Il1bnegS100a9neg" "31" "42"

############ 32/42 template_Metadatafiltered_Il1bnegS100a9neg.R ########
run_workbook_template "template_Metadatafiltered_Il1bnegS100a9neg" "32" "42"

############ 33/42 template_ReclusterSO_Il1bnegS100a9neg.R ########
run_workbook_template "template_ReclusterSO_Il1bnegS100a9neg" "33" "42"

############ 34/42 template_SampleNames_Il1bnegS100a9neg.R ########
run_workbook_template "template_SampleNames_Il1bnegS100a9neg" "34" "42"

############ 35/42 template_MetadataTable_recluster_Il1bnegS100a9neg.R ########
run_workbook_template "template_MetadataTable_recluster_Il1bnegS100a9neg" "35" "42"

############ 36/42 template_DEGMarkers_Il1bnegS100a9neg.R ########
run_workbook_template "template_DEGMarkers_Il1bnegS100a9neg" "36" "42"

############ 37/42 template_ENTN803aPD1vsENTN803.R ########
run_workbook_template "template_ENTN803aPD1vsENTN803" "37" "42"

############ 38/42 template_FIGURE_4D.R ########
run_workbook_template "template_FIGURE_4D" "38" "42"

############ 39/42 template_GSEA_Preranked1.R ########
run_workbook_template "template_GSEA_Preranked1" "39" "42"

############ 40/42 template_L2PSingle_2.R ########
run_workbook_template "template_L2PSingle_2" "40" "42"

############ 41/42 template_GSEA_Filtered1.R ########
run_workbook_template "template_GSEA_Filtered1" "41" "42"

############ 42/42 template_FIGURE_4E.R ########
run_workbook_template "template_FIGURE_4E" "42" "42"


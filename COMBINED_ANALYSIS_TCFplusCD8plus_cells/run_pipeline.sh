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
############ 1/52 template_PostIt.R ########
run_workbook_template "template_PostIt" "1" "52"

############ 2/52 template_PostIt_1.R ########
run_workbook_template "template_PostIt_1" "2" "52"

############ 3/52 template_PostIt_2.R ########
run_workbook_template "template_PostIt_2" "3" "52"

############ 4/52 template_PostIt_CellTypes.R ########
run_workbook_template "template_PostIt_CellTypes" "4" "52"

############ 5/52 template_PostIt_CombNorm.R ########
run_workbook_template "template_PostIt_CombNorm" "5" "52"

############ 6/52 template_PostIt_DualLabel.R ########
run_workbook_template "template_PostIt_DualLabel" "6" "52"

############ 7/52 template_PostIt_Filtered.R ########
run_workbook_template "template_PostIt_Filtered" "7" "52"

############ 8/52 template_PostIt_ProcessedInput.R ########
run_workbook_template "template_PostIt_ProcessedInput" "8" "52"

############ 9/52 template_PostIt_QCFiltered.R ########
run_workbook_template "template_PostIt_QCFiltered" "9" "52"

############ 10/52 template_PostIt_Reclustered.R ########
run_workbook_template "template_PostIt_Reclustered" "10" "52"

############ 11/52 template_PostIt_YourDataHere.R ########
run_workbook_template "template_PostIt_YourDataHere" "11" "52"

############ 12/52 template_ProcessedInputSO.R ########
run_workbook_template "template_ProcessedInputSO" "12" "52"

############ 13/52 template_QCFilteredSO.R ########
run_workbook_template "template_QCFilteredSO" "13" "52"

############ 14/52 template_CombNormSO.R ########
run_workbook_template "template_CombNormSO" "14" "52"

############ 15/52 template_MetadataTable_CombNorm.R ########
run_workbook_template "template_MetadataTable_CombNorm" "15" "52"

############ 16/52 template_CellTypesSO.R ########
run_workbook_template "template_CellTypesSO" "16" "52"

############ 17/52 template_MetadataTable_CellTypes.R ########
run_workbook_template "template_MetadataTable_CellTypes" "17" "52"

############ 18/52 template_SampleNames_CellTypes.R ########
run_workbook_template "template_SampleNames_CellTypes" "18" "52"

############ 19/52 template_DualLabelSO_CD3CD8_Tcells.R ########
run_workbook_template "template_DualLabelSO_CD3CD8_Tcells" "19" "52"

############ 20/52 template_FilteredSO.R ########
run_workbook_template "template_FilteredSO" "20" "52"

############ 21/52 template_MDT_DL_CD3CD8_Tcells.R ########
run_workbook_template "template_MDT_DL_CD3CD8_Tcells" "21" "52"

############ 22/52 template_MetadataTable_Filtered_9_10.R ########
run_workbook_template "template_MetadataTable_Filtered_9_10" "22" "52"

############ 23/52 template_ReclusteredSO.R ########
run_workbook_template "template_ReclusteredSO" "23" "52"

############ 24/52 template_SampleNames_Reclustered.R ########
run_workbook_template "template_SampleNames_Reclustered" "24" "52"

############ 25/52 template_FSO_CD3CD8_Tcells.R ########
run_workbook_template "template_FSO_CD3CD8_Tcells" "25" "52"

############ 26/52 template_MDT_FSO_CD3CD8_Tcells.R ########
run_workbook_template "template_MDT_FSO_CD3CD8_Tcells" "26" "52"

############ 27/52 template_MetadataTable_Reclustered.R ########
run_workbook_template "template_MetadataTable_Reclustered" "27" "52"

############ 28/52 template_Reclustered_CD3CD8_Tcells.R ########
run_workbook_template "template_Reclustered_CD3CD8_Tcells" "28" "52"

############ 29/52 template_TCF7pos_Tim3neg.R ########
run_workbook_template "template_TCF7pos_Tim3neg" "29" "52"

############ 30/52 template_DEGMarkers.R ########
run_workbook_template "template_DEGMarkers" "30" "52"

############ 31/52 template_MDT_TCF7pos_Tim3neg.R ########
run_workbook_template "template_MDT_TCF7pos_Tim3neg" "31" "52"

############ 32/52 template_RMDT_CD3CD8_Tcells.R ########
run_workbook_template "template_RMDT_CD3CD8_Tcells" "32" "52"

############ 33/52 template_RSN_CD3CD8_Tcells.R ########
run_workbook_template "template_RSN_CD3CD8_Tcells" "33" "52"

############ 34/52 template_DEGMarkers_1.R ########
run_workbook_template "template_DEGMarkers_1" "34" "52"

############ 35/52 template_FSO_TCF7pos_Tim3neg.R ########
run_workbook_template "template_FSO_TCF7pos_Tim3neg" "35" "52"

############ 36/52 template_DLSO_Cd44posGzmbneg.R ########
run_workbook_template "template_DLSO_Cd44posGzmbneg" "36" "52"

############ 37/52 template_FMDT_TCF7pos_Tim3neg.R ########
run_workbook_template "template_FMDT_TCF7pos_Tim3neg" "37" "52"

############ 38/52 template_RSO_TCF7pos_Tim3neg.R ########
run_workbook_template "template_RSO_TCF7pos_Tim3neg" "38" "52"

############ 39/52 template_DLMDT_Cd44posGzmbneg.R ########
run_workbook_template "template_DLMDT_Cd44posGzmbneg" "39" "52"

############ 40/52 template_RMDT_TCF7pos_Tim3neg.R ########
run_workbook_template "template_RMDT_TCF7pos_Tim3neg" "40" "52"

############ 41/52 template_RSN_TCF7pos_Tim3neg.R ########
run_workbook_template "template_RSN_TCF7pos_Tim3neg" "41" "52"

############ 42/52 template_TCF7pos_Tim3neg_DEG.R ########
run_workbook_template "template_TCF7pos_Tim3neg_DEG" "42" "52"

############ 43/52 template_FIGURE_3N_TripletvsENTN803.R ########
run_workbook_template "template_FIGURE_3N_TripletvsENTN803" "43" "52"

############ 44/52 template_FIGURE_3N_TripletvsENTaPD1.R ########
run_workbook_template "template_FIGURE_3N_TripletvsENTaPD1" "44" "52"

############ 45/52 template_FIGURE_3N_TripletvsPBS.R ########
run_workbook_template "template_FIGURE_3N_TripletvsPBS" "45" "52"

############ 46/52 template_FSO_Cd44posGzmbneg.R ########
run_workbook_template "template_FSO_Cd44posGzmbneg" "46" "52"

############ 47/52 template_FMDT_Cd44posGzmbneg.R ########
run_workbook_template "template_FMDT_Cd44posGzmbneg" "47" "52"

############ 48/52 template_RSO_Cd44posGzmbneg.R ########
run_workbook_template "template_RSO_Cd44posGzmbneg" "48" "52"

############ 49/52 template_RMDT_Cd44posGzmbneg.R ########
run_workbook_template "template_RMDT_Cd44posGzmbneg" "49" "52"

############ 50/52 template_RSN_Cd44posGzmbneg.R ########
run_workbook_template "template_RSN_Cd44posGzmbneg" "50" "52"

############ 51/52 template_Cd44posGzmbneg_DEG.R ########
run_workbook_template "template_Cd44posGzmbneg_DEG" "51" "52"

############ 52/52 template_FIGURE_3J.R ########
run_workbook_template "template_FIGURE_3J" "52" "52"


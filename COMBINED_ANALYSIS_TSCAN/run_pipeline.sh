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
############ 1/46 template_PostIt.R ########
run_workbook_template "template_PostIt" "1" "46"

############ 2/46 template_PostIt_CellTypes.R ########
run_workbook_template "template_PostIt_CellTypes" "2" "46"

############ 3/46 template_PostIt_CombNorm.R ########
run_workbook_template "template_PostIt_CombNorm" "3" "46"

############ 4/46 template_PostIt_DualLabel.R ########
run_workbook_template "template_PostIt_DualLabel" "4" "46"

############ 5/46 template_PostIt_Filtered.R ########
run_workbook_template "template_PostIt_Filtered" "5" "46"

############ 6/46 template_PostIt_ProcessedInput.R ########
run_workbook_template "template_PostIt_ProcessedInput" "6" "46"

############ 7/46 template_PostIt_QCFiltered.R ########
run_workbook_template "template_PostIt_QCFiltered" "7" "46"

############ 8/46 template_PostIt_YourDataHere.R ########
run_workbook_template "template_PostIt_YourDataHere" "8" "46"

############ 9/46 template_PostIt_main.R ########
run_workbook_template "template_PostIt_main" "9" "46"

############ 10/46 template_ProcessedInputSO.R ########
run_workbook_template "template_ProcessedInputSO" "10" "46"

############ 11/46 template_QCFilteredSO.R ########
run_workbook_template "template_QCFilteredSO" "11" "46"

############ 12/46 template_CombNormSO.R ########
run_workbook_template "template_CombNormSO" "12" "46"

############ 13/46 template_MetadataTable_CombNorm.R ########
run_workbook_template "template_MetadataTable_CombNorm" "13" "46"

############ 14/46 template_CellTypesSO.R ########
run_workbook_template "template_CellTypesSO" "14" "46"

############ 15/46 template_MetadataTable_CellTypes.R ########
run_workbook_template "template_MetadataTable_CellTypes" "15" "46"

############ 16/46 template_SampleNames_CellTypes.R ########
run_workbook_template "template_SampleNames_CellTypes" "16" "46"

############ 17/46 template_DualLabelSO_CD3CD8_Tcells.R ########
run_workbook_template "template_DualLabelSO_CD3CD8_Tcells" "17" "46"

############ 18/46 template_FilteredSO.R ########
run_workbook_template "template_FilteredSO" "18" "46"

############ 19/46 template_MDT_DL_CD3CD8_Tcells.R ########
run_workbook_template "template_MDT_DL_CD3CD8_Tcells" "19" "46"

############ 20/46 template_FSO_CD3CD8_Tcells.R ########
run_workbook_template "template_FSO_CD3CD8_Tcells" "20" "46"

############ 21/46 template_MDT_FSO_CD3CD8_Tcells.R ########
run_workbook_template "template_MDT_FSO_CD3CD8_Tcells" "21" "46"

############ 22/46 template_ExcludeNK_CD8Tcells.R ########
run_workbook_template "template_ExcludeNK_CD8Tcells" "22" "46"

############ 23/46 template_MDT_ExcludeNK_CD8Tcells.R ########
run_workbook_template "template_MDT_ExcludeNK_CD8Tcells" "23" "46"

############ 24/46 template_FSO_ExcludeNK_CD8Tcells.R ########
run_workbook_template "template_FSO_ExcludeNK_CD8Tcells" "24" "46"

############ 25/46 template_FMDT_ExcludeNK_CD8Tcells.R ########
run_workbook_template "template_FMDT_ExcludeNK_CD8Tcells" "25" "46"

############ 26/46 template_RSO_ExcludeNK_CD8Tcells.R ########
run_workbook_template "template_RSO_ExcludeNK_CD8Tcells" "26" "46"

############ 27/46 template_ColorByMeta.R ########
run_workbook_template "template_ColorByMeta" "27" "46"

############ 28/46 template_RMDT_ExcludeNK_CD8Tcells.R ########
run_workbook_template "template_RMDT_ExcludeNK_CD8Tcells" "28" "46"

############ 29/46 template_RSN_ExcludeNK_CD8Tcells.R ########
run_workbook_template "template_RSN_ExcludeNK_CD8Tcells" "29" "46"

############ 30/46 template_SubsettedSO.R ########
run_workbook_template "template_SubsettedSO" "30" "46"

############ 31/46 template_TrajectoryAnalysis.R ########
run_workbook_template "template_TrajectoryAnalysis" "31" "46"

############ 32/46 template_ColorByGene_1.R ########
run_workbook_template "template_ColorByGene_1" "32" "46"

############ 33/46 template_ColorByGene_4.R ########
run_workbook_template "template_ColorByGene_4" "33" "46"

############ 34/46 template_ColorByMeta_1.R ########
run_workbook_template "template_ColorByMeta_1" "34" "46"

############ 35/46 template_ColorByMeta_2.R ########
run_workbook_template "template_ColorByMeta_2" "35" "46"

############ 36/46 template_ColorByMetadata_1.R ########
run_workbook_template "template_ColorByMetadata_1" "36" "46"

############ 37/46 template_ColorByMetadata_3.R ########
run_workbook_template "template_ColorByMetadata_3" "37" "46"

############ 38/46 template_DEGMarkers.R ########
run_workbook_template "template_DEGMarkers" "38" "46"

############ 39/46 template_MetadataTable_7.R ########
run_workbook_template "template_MetadataTable_7" "39" "46"

############ 40/46 template_ReclusterSO_1.R ########
run_workbook_template "template_ReclusterSO_1" "40" "46"

############ 41/46 template_SampleNames_2.R ########
run_workbook_template "template_SampleNames_2" "41" "46"

############ 42/46 template_MetadataTable_8.R ########
run_workbook_template "template_MetadataTable_8" "42" "46"

############ 43/46 template_Traj_2.R ########
run_workbook_template "template_Traj_2" "43" "46"

############ 44/46 template_ColorByMeta_3.R ########
run_workbook_template "template_ColorByMeta_3" "44" "46"

############ 45/46 template_ColorByMeta_4.R ########
run_workbook_template "template_ColorByMeta_4" "45" "46"

############ 46/46 template_ColorByMeta_5.R ########
run_workbook_template "template_ColorByMeta_5" "46" "46"


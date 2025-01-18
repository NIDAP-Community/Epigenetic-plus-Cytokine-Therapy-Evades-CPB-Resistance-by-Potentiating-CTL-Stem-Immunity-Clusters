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
############ 1/68 template_PostIt_CellTypes.R ########
run_workbook_template "template_PostIt_CellTypes" "1" "68"

############ 2/68 template_PostIt_CombNorm.R ########
run_workbook_template "template_PostIt_CombNorm" "2" "68"

############ 3/68 template_PostIt_DualLabel.R ########
run_workbook_template "template_PostIt_DualLabel" "3" "68"

############ 4/68 template_PostIt_Filtered.R ########
run_workbook_template "template_PostIt_Filtered" "4" "68"

############ 5/68 template_PostIt_ProcessedInput.R ########
run_workbook_template "template_PostIt_ProcessedInput" "5" "68"

############ 6/68 template_PostIt_QCFiltered.R ########
run_workbook_template "template_PostIt_QCFiltered" "6" "68"

############ 7/68 template_PostIt_Reclustered.R ########
run_workbook_template "template_PostIt_Reclustered" "7" "68"

############ 8/68 template_PostIt_YourDataHere.R ########
run_workbook_template "template_PostIt_YourDataHere" "8" "68"

############ 9/68 template_ProcessedInputSO.R ########
run_workbook_template "template_ProcessedInputSO" "9" "68"

############ 10/68 template_QCFilteredSO.R ########
run_workbook_template "template_QCFilteredSO" "10" "68"

############ 11/68 template_CombNormSO.R ########
run_workbook_template "template_CombNormSO" "11" "68"

############ 12/68 template_MetadataTable_CombNorm.R ########
run_workbook_template "template_MetadataTable_CombNorm" "12" "68"

############ 13/68 template_CellTypesSO.R ########
run_workbook_template "template_CellTypesSO" "13" "68"

############ 14/68 template_MetadataTable_CellTypes.R ########
run_workbook_template "template_MetadataTable_CellTypes" "14" "68"

############ 15/68 template_SampleNames_CellTypes.R ########
run_workbook_template "template_SampleNames_CellTypes" "15" "68"

############ 16/68 template_DualLabelSO_CD8CD3.R ########
run_workbook_template "template_DualLabelSO_CD8CD3" "16" "68"

############ 17/68 template_FilteredSO.R ########
run_workbook_template "template_FilteredSO" "17" "68"

############ 18/68 template_MetadataTableCD8CD3.R ########
run_workbook_template "template_MetadataTableCD8CD3" "18" "68"

############ 19/68 template_MetadataTable_Filtered_9_10.R ########
run_workbook_template "template_MetadataTable_Filtered_9_10" "19" "68"

############ 20/68 template_ReclusteredSO_CD8Tcells.R ########
run_workbook_template "template_ReclusteredSO_CD8Tcells" "20" "68"

############ 21/68 template_SampleNames_Reclustered.R ########
run_workbook_template "template_SampleNames_Reclustered" "21" "68"

############ 22/68 template_FilteredSO_CD8CD3.R ########
run_workbook_template "template_FilteredSO_CD8CD3" "22" "68"

############ 23/68 template_MetadataTable_CD8cd3.R ########
run_workbook_template "template_MetadataTable_CD8cd3" "23" "68"

############ 24/68 template_MetadataTable_Reclustered.R ########
run_workbook_template "template_MetadataTable_Reclustered" "24" "68"

############ 25/68 template_ReclusterSO_CD8cd3.R ########
run_workbook_template "template_ReclusterSO_CD8cd3" "25" "68"

############ 26/68 template_SampleNames_Reclustered_CD8CD3.R ########
run_workbook_template "template_SampleNames_Reclustered_CD8CD3" "26" "68"

############ 27/68 template_SuppFIGURE_3A_ENT.R ########
run_workbook_template "template_SuppFIGURE_3A_ENT" "27" "68"

############ 28/68 template_SuppFIGURE_3A_ENTN803.R ########
run_workbook_template "template_SuppFIGURE_3A_ENTN803" "28" "68"

############ 29/68 template_SuppFIGURE_3A_ENTN803aPD1.R ########
run_workbook_template "template_SuppFIGURE_3A_ENTN803aPD1" "29" "68"

############ 30/68 template_SuppFIGURE_3A_ENTaPD1.R ########
run_workbook_template "template_SuppFIGURE_3A_ENTaPD1" "30" "68"

############ 31/68 template_SuppFIGURE_3A_N803.R ########
run_workbook_template "template_SuppFIGURE_3A_N803" "31" "68"

############ 32/68 template_SuppFIGURE_3A_N803aPD1.R ########
run_workbook_template "template_SuppFIGURE_3A_N803aPD1" "32" "68"

############ 33/68 template_SuppFIGURE_3A_PBS.R ########
run_workbook_template "template_SuppFIGURE_3A_PBS" "33" "68"

############ 34/68 template_SuppFIGURE_3A_aPD1.R ########
run_workbook_template "template_SuppFIGURE_3A_aPD1" "34" "68"

############ 35/68 template_SuppFIGURE_4A.R ########
run_workbook_template "template_SuppFIGURE_4A" "35" "68"

############ 36/68 template_DEGMarkers.R ########
run_workbook_template "template_DEGMarkers" "36" "68"

############ 37/68 template_DualLabelSO_CD8T_Gzmb_Mki67neg_copied.R ########
run_workbook_template "template_DualLabelSO_CD8T_Gzmb_Mki67neg_copied" "37" "68"

############ 38/68 template_MetadataTable_CD8T_Gzmb_Mki67neg_copied.R ########
run_workbook_template "template_MetadataTable_CD8T_Gzmb_Mki67neg_copied" "38" "68"

############ 39/68 template_MetadataTable_Reclustered_CD8CD3.R ########
run_workbook_template "template_MetadataTable_Reclustered_CD8CD3" "39" "68"

############ 40/68 template_DEGMarkers_Reclustered_CD8CD3.R ########
run_workbook_template "template_DEGMarkers_Reclustered_CD8CD3" "40" "68"

############ 41/68 template_FIGURE_2I.R ########
run_workbook_template "template_FIGURE_2I" "41" "68"

############ 42/68 template_FIGURE_2J.R ########
run_workbook_template "template_FIGURE_2J" "42" "68"

############ 43/68 template_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg.R ########
run_workbook_template "template_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg" "43" "68"

############ 44/68 template_MetadataTable_CD8T_Gzmb_Mki67_copied.R ########
run_workbook_template "template_MetadataTable_CD8T_Gzmb_Mki67_copied" "44" "68"

############ 45/68 template_ReclusterSO_CD8_GzmbposKi67neg.R ########
run_workbook_template "template_ReclusterSO_CD8_GzmbposKi67neg" "45" "68"

############ 46/68 template_SampleNames_CD8_GzmbposKi67neg.R ########
run_workbook_template "template_SampleNames_CD8_GzmbposKi67neg" "46" "68"

############ 47/68 template_SuppFIGURE_2F_ENT_vs_PBS.R ########
run_workbook_template "template_SuppFIGURE_2F_ENT_vs_PBS" "47" "68"

############ 48/68 template_SuppFIGURE_2F_N803_vs_PBS.R ########
run_workbook_template "template_SuppFIGURE_2F_N803_vs_PBS" "48" "68"

############ 49/68 template_SuppFIGURE_2F_Triplet_vs_ENTaN803.R ########
run_workbook_template "template_SuppFIGURE_2F_Triplet_vs_ENTaN803" "49" "68"

############ 50/68 template_SuppFIGURE_2F_Triplet_vs_ENTaPD1.R ########
run_workbook_template "template_SuppFIGURE_2F_Triplet_vs_ENTaPD1" "50" "68"

############ 51/68 template_SuppFIGURE_2F_Triplet_vs_N803aPD1.R ########
run_workbook_template "template_SuppFIGURE_2F_Triplet_vs_N803aPD1" "51" "68"

############ 52/68 template_SuppFIGURE_2F_aPD1_vs_PBS.R ########
run_workbook_template "template_SuppFIGURE_2F_aPD1_vs_PBS" "52" "68"

############ 53/68 template_SuppFIGURE_2G_Triplet_vs_ENTN803_FINAL_ADDTOGEO.R ########
run_workbook_template "template_SuppFIGURE_2G_Triplet_vs_ENTN803_FINAL_ADDTOGEO" "53" "68"

############ 54/68 template_SuppFIGURE_2G_Triplet_vs_ENTN803_OLD_DONOTADD.R ########
run_workbook_template "template_SuppFIGURE_2G_Triplet_vs_ENTN803_OLD_DONOTADD" "54" "68"

############ 55/68 template_SuppFIGURE_2G_Triplet_vs_ENTaPD1.R ########
run_workbook_template "template_SuppFIGURE_2G_Triplet_vs_ENTaPD1" "55" "68"

############ 56/68 template_MetadataTable_CD8GzmbposKi67neg.R ########
run_workbook_template "template_MetadataTable_CD8GzmbposKi67neg" "56" "68"

############ 57/68 template_DEGMarkers_CD8_GzmbposKi67neg.R ########
run_workbook_template "template_DEGMarkers_CD8_GzmbposKi67neg" "57" "68"

############ 58/68 template_ENTN803aPD1_vs_PBS.R ########
run_workbook_template "template_ENTN803aPD1_vs_PBS" "58" "68"

############ 59/68 template_GSEA_Preranked.R ########
run_workbook_template "template_GSEA_Preranked" "59" "68"

############ 60/68 template_SuppFIGURE_3B.R ########
run_workbook_template "template_SuppFIGURE_3B" "60" "68"

############ 61/68 template_SuppFIGURE_3C_ENTvsPBS.R ########
run_workbook_template "template_SuppFIGURE_3C_ENTvsPBS" "61" "68"

############ 62/68 template_SuppFIGURE_3C_N803vsPBS.R ########
run_workbook_template "template_SuppFIGURE_3C_N803vsPBS" "62" "68"

############ 63/68 template_SuppFIGURE_3C_TripletvsENTN803.R ########
run_workbook_template "template_SuppFIGURE_3C_TripletvsENTN803" "63" "68"

############ 64/68 template_SuppFIGURE_3C_TripletvsENTaPD1.R ########
run_workbook_template "template_SuppFIGURE_3C_TripletvsENTaPD1" "64" "68"

############ 65/68 template_SuppFIGURE_3C_TripletvsN803aPD1.R ########
run_workbook_template "template_SuppFIGURE_3C_TripletvsN803aPD1" "65" "68"

############ 66/68 template_SuppFIGURE_3C_TripletvsPBS.R ########
run_workbook_template "template_SuppFIGURE_3C_TripletvsPBS" "66" "68"

############ 67/68 template_SuppFIGURE_3C_aPD1vsPBS.R ########
run_workbook_template "template_SuppFIGURE_3C_aPD1vsPBS" "67" "68"

############ 68/68 template_GSEAFiltered_CD8GzmbposKi67neg.R ########
run_workbook_template "template_GSEAFiltered_CD8GzmbposKi67neg" "68" "68"


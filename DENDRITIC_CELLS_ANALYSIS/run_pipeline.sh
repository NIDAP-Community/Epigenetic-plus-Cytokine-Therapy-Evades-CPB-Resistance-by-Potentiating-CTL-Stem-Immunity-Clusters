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
############ 1/37 template_PostIt.R ########
run_workbook_template "template_PostIt" "1" "37"

############ 2/37 template_PostIt_1.R ########
run_workbook_template "template_PostIt_1" "2" "37"

############ 3/37 template_PostIt_CombNorm.R ########
run_workbook_template "template_PostIt_CombNorm" "3" "37"

############ 4/37 template_PostIt_DualLabel.R ########
run_workbook_template "template_PostIt_DualLabel" "4" "37"

############ 5/37 template_PostIt_ProcessedInput.R ########
run_workbook_template "template_PostIt_ProcessedInput" "5" "37"

############ 6/37 template_PostIt_QCFiltered.R ########
run_workbook_template "template_PostIt_QCFiltered" "6" "37"

############ 7/37 template_PostIt_YourDataHere.R ########
run_workbook_template "template_PostIt_YourDataHere" "7" "37"

############ 8/37 template_ProcessedInputSO.R ########
run_workbook_template "template_ProcessedInputSO" "8" "37"

############ 9/37 template_QCFilteredSO.R ########
run_workbook_template "template_QCFilteredSO" "9" "37"

############ 10/37 template_CombNormSO.R ########
run_workbook_template "template_CombNormSO" "10" "37"

############ 11/37 template_MetadataTable_CombNorm.R ########
run_workbook_template "template_MetadataTable_CombNorm" "11" "37"

############ 12/37 template_CellTypesSO.R ########
run_workbook_template "template_CellTypesSO" "12" "37"

############ 13/37 template_MetadataTable_CellTypes.R ########
run_workbook_template "template_MetadataTable_CellTypes" "13" "37"

############ 14/37 template_SampleNames_CellTypes.R ########
run_workbook_template "template_SampleNames_CellTypes" "14" "37"

############ 15/37 template_DualLabelSO_CD11cNcr1.R ########
run_workbook_template "template_DualLabelSO_CD11cNcr1" "15" "37"

############ 16/37 template_FilteredSO.R ########
run_workbook_template "template_FilteredSO" "16" "37"

############ 17/37 template_MetadataTable_DualLabelCD11cNcr1.R ########
run_workbook_template "template_MetadataTable_DualLabelCD11cNcr1" "17" "37"

############ 18/37 template_Filtered_SO_DualLabel_CD11cNcr1.R ########
run_workbook_template "template_Filtered_SO_DualLabel_CD11cNcr1" "18" "37"

############ 19/37 template_Metadata_Table_DualLabel_CD11cNcr1.R ########
run_workbook_template "template_Metadata_Table_DualLabel_CD11cNcr1" "19" "37"

############ 20/37 template_ReclusterSO_DualLabel_CD11cNcr1.R ########
run_workbook_template "template_ReclusterSO_DualLabel_CD11cNcr1" "20" "37"

############ 21/37 template_SampleNames_DualLabel_CD11cNcr1.R ########
run_workbook_template "template_SampleNames_DualLabel_CD11cNcr1" "21" "37"

############ 22/37 template_Clec9aItgaeDP.R ########
run_workbook_template "template_Clec9aItgaeDP" "22" "37"

############ 23/37 template_MetadataTable_reclustered_DualLabel_CD11cNcr1.R ########
run_workbook_template "template_MetadataTable_reclustered_DualLabel_CD11cNcr1" "23" "37"

############ 24/37 template_Metadata_Clec9aItgaeDP.R ########
run_workbook_template "template_Metadata_Clec9aItgaeDP" "24" "37"

############ 25/37 template_FilteredSO_Clec9aItgaeDP.R ########
run_workbook_template "template_FilteredSO_Clec9aItgaeDP" "25" "37"

############ 26/37 template_Filteredmetadata_Clec9aItgaeDP.R ########
run_workbook_template "template_Filteredmetadata_Clec9aItgaeDP" "26" "37"

############ 27/37 template_ReclusterSO_Clec9aItgaeDP.R ########
run_workbook_template "template_ReclusterSO_Clec9aItgaeDP" "27" "37"

############ 28/37 template_Reclusteredmetadata_Clec9aItgaeDP.R ########
run_workbook_template "template_Reclusteredmetadata_Clec9aItgaeDP" "28" "37"

############ 29/37 template_Samplenames_Clec9aItgaeDP.R ########
run_workbook_template "template_Samplenames_Clec9aItgaeDP" "29" "37"

############ 30/37 template_DEGmarkers_Clec9aItgaeDP.R ########
run_workbook_template "template_DEGmarkers_Clec9aItgaeDP" "30" "37"

############ 31/37 template_FIGURE_5A.R ########
run_workbook_template "template_FIGURE_5A" "31" "37"

############ 32/37 template_FIGURE_5C_Triplet_vs_PBS.R ########
run_workbook_template "template_FIGURE_5C_Triplet_vs_PBS" "32" "37"

############ 33/37 template_FIGURE_5D_ENTN803aPD1_vs_PBS.R ########
run_workbook_template "template_FIGURE_5D_ENTN803aPD1_vs_PBS" "33" "37"

############ 34/37 template_GSEA_Preranked.R ########
run_workbook_template "template_GSEA_Preranked" "34" "37"

############ 35/37 template_GSEA_Filtered.R ########
run_workbook_template "template_GSEA_Filtered" "35" "37"

############ 36/37 template_FIGURE_5Bpart1.R ########
run_workbook_template "template_FIGURE_5Bpart1" "36" "37"

############ 37/37 template_FIGURE_5Bpart2.R ########
run_workbook_template "template_FIGURE_5Bpart2" "37" "37"


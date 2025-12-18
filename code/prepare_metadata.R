# prepare_metadata.R

library(dplyr)
library(pals) 
library(data.table)
library(readxl)

# Directories
project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
processed_dir <- file.path(data_dir, "analyzed")

# Load Data (From previous step)
message("   Loading Metadata...")
targets <- fread(file.path(processed_dir, "GSE140686_Final_Targets.csv"))
meta_file <- file.path(data_dir, "raw/GSE140686/41467_2020_20603_MOESM4_ESM.xlsx")
meta_df <- read_excel(meta_file)

message("   Generating Methylation Class Acronyms (Methyl_Dx)...")

methyl_map <- tribble(
  ~Methylation_Class, ~Methyl_Dx,
  "rhabdomyosarcoma (alveolar)", "RMS-A",
  "rhabdomyosarcoma (embryonal)", "RMS-E",
  "synovial sarcoma", "SS",
  "myxoid liposarcoma", "MLPS",
  "small blue round cell tumour with BCOR alteration", "BCOR",
  "small blue round cell tumour with CIC alteration", "CIC",
  "Ewing´s sarcoma", "EWS",
  "fibrous dysplasia", "FD",
  "chordoma", "CHORD",
  "osteosarcoma (high grade)", "OS-HG",
  "undifferentiated sarcoma", "US",
  "chondrosarcoma (IDH group A)", "CSA-IDH-A",
  "chondrosarcoma (group B)", "CSA-B",
  "chondrosarcoma (mesenchymal)", "CSA-M",
  "chondrosarcoma (group A)", "CSA-A",
  "extraskeletal myxoid chondrosarcoma", "EMC",
  "desmoplastic small round cell tumour", "DSRCT",
  "alveolar soft part sarcoma", "ASPS",
  "sarcoma (RMS-like)", "SARC-RMS",
  "solitary fibrous tumour", "SFT",
  "sarcoma (MPNST-like)", "SARC-MPNST",
  "angiosarcoma", "AS",
  "gastrointestinal stromal tumour", "GIST",
  "leiomyosarcoma", "LMS",
  "control (reactive tissue)", "CTRL-TISSUE",
  "chordoma (dedifferentiated)", "CHORD-DD",
  "lipoma", "LPO",
  "dermatofibrosarcoma protuberans", "DFSP",
  "giant cell tumour of bone", "GCTB",
  "malignant peripheral nerve sheath tumour", "MPNST",
  "osteoblastoma", "OBL",
  "desmoid-type fibromatosis", "DES",
  "low-grade fibromyxoid sarcoma", "LGFMS",
  "neurofibroma", "NF",
  "neurofibroma (plexiform)", "NF-P",
  "schwannoma", "SCHW",
  "infantile fibrosarcoma", "IFS",
  "rhabdomyosarcoma (MYOD1)", "RMS-MYOD1",
  "Kaposi sarcoma", "KS",
  "melanoma (cutaneous)", "MEL",
  "sclerosing epithelioid fibrosarcoma", "SEF",
  "well- / dedifferentiated liposarcoma", "WDLPS/DDLPS",
  "control (muscle tissue)", "CTRL-MUSCLE",
  "malignant rhabdoid tumour", "MRT",
  "chondrosarcoma (clear cell)", "CSA-CC",
  "chondroblastoma", "CB",
  "chondrosarcoma (IDH group B)", "CSA-IDH-B",
  "squamous cell carcinoma (cutaneous)", "SCC",
  "endometrial stromal sarcoma (low grade)", "ESS-LG",
  "atypical fibroxanthoma / pleomorphic dermal sarcoma", "AFX/PDS",
  "clear cell sarcoma of the kidney", "CCS-K",
  "epithelioid sarcoma", "ES",
  "myositis ossificans", "MOS",
  "epithelioid haemangioendothelioma", "EHE",
  "inflammatory myofibroblastic tumour", "IMT",
  "nodular fasciitis", "NOD",
  "angiomatoid fibrous histiocytoma", "AFH",
  "myositis proliferans", "MYP",
  "ossifying fibromyxoid tumour", "OFMT",
  "Langerhans cell histiocytosis", "LCH",
  "leiomyoma", "LMO",
  "clear cell sarcoma of soft parts", "CCS",
  "angioleiomyoma / myopericytoma", "ALMO/MPC",
  "endometrial stromal sarcoma (high grade)", "ESS-HG",
  "control (blood)", "CTRL-BLOOD"
)
    
# Merge with targets and had # to colors
temp_meta <- targets %>%
    left_join(methyl_map, by = "Methylation_Class") %>% 
    mutate(
      Color = paste0("#", Color), 
      ID = paste(geo_accession, IDAT, sep = "_")) %>%
    rename(Color_Methyl = Color)

message("   Adding diagnosis acronyms and color palette...")

dx_map <- tribble(
  ~Diagnosis, ~Dx, ~WHO_differentiation,
  'Angioleiomyoma', 'ALMO', 'Perivascular tumors',
  'Leiomyoma', 'LMO', 'Smooth muscle tumors',
  'Leiomyosarcoma', 'LMS', 'Smooth muscle tumors',
  'Myopericytoma', 'MPC', 'Perivascular tumors',
  'Angiomatoid fibrous histiocytoma', 'AFH', 'Tumors of uncertain differentiation',
  'Atypical fibroxanthoma', 'AFX', 'Tumors of uncertain differentiation',
  'Rhabdomyosarcoma (alveolar)', 'RMS-A', 'Skeletal muscle tumors',
  'Rhabdomyosarcoma (embryonal)', 'RMS-E', 'Skeletal muscle tumors',
  'Embryonal rhabdomyosarcoma', 'RMS-E', 'Skeletal muscle tumors',
  'Rhabdomyosarcoma (spindle cell)', 'RMS-S', 'Skeletal muscle tumors',
  'Rhabdomyosarcoma (sclerosing)', 'RMS-S', 'Skeletal muscle tumors',
  'Malignant rhabdoid tumour', 'MRT', 'Tumors of uncertain differentiation',
  'Synovial sarcoma', 'SS', 'Tumors of uncertain differentiation',
  'Myxoid liposarcoma', 'MLPS', 'Adipocytic tumors',
  'Small blue round cell tumour with BCOR alteration', 'BCOR', 'Undifferentiated small round cell tumors',
  'Small blue round cell tumour with CIC alteration', 'CIC', 'Undifferentiated small round cell tumors',
  'Ewing´s sarcoma', 'EWS', 'Undifferentiated small round cell tumors',
  'Fibrous dysplasia', 'FD', 'Other mesenchymal tumors of bone',
  'Chordoma', 'CHORD', 'Notochordal tumors',
  'Osteosarcoma (high-grade)', 'OS-HG', 'Osteogenic tumors',
  'Undifferentiated pleomorphic sarcoma', 'UPS', 'Undifferentiated sarcomas',
  'Pleomorphic dermal sarcoma', 'PDS', 'Undifferentiated sarcomas',
  'Undifferentiated epithelioid sarcoma not otherwise specified', 'UES-NOS', 'Undifferentiated sarcomas',
  'Chondrosarcoma', 'CSA', 'Chondrogenic tumors',
  'Chondrosarcoma (mesenchymal)', 'CSA-M', 'Chondrogenic tumors',
  'Mesenchymal chondrosarcoma', 'CSA-M', 'Chondrogenic tumors',
  'Chondrosarcoma (clear cell)', 'CSA-CC', 'Chondrogenic tumors',
  'Chondroblastoma', 'CB', 'Chondrogenic tumors',
  'Extraskeletal myxoid chondrosarcoma', 'EMC', 'Tumors of uncertain differentiation',
  'Clear cell sarcoma (kidney)', 'CCS-K', 'Tumors of uncertain differentiation',
  'Clear cell sarcoma (soft tissue)', 'CCS', 'Tumors of uncertain differentiation',
  'Desmoplastic small round cell tumour', 'DSRCT', 'Tumors of uncertain differentiation',
  'Alveolar soft part sarcoma', 'ASPS', 'Tumors of uncertain differentiation',
  'Primitive neuroectodermal tumour', 'EWS', 'Undifferentiated small round cell tumors',
  'Solitary fibrous tumour', 'SFT', 'Fibroblastic/myofibroblastic tumors',
  'Malignant peripheral nerve sheath tumour', 'MPNST', 'Peripheral nerve sheath tumors',
  'Angiosarcoma', 'AS', 'Vascular tumors',
  'Gastroinstestinal stromal tumour', 'GIST', 'Gastrointestinal stromal tumor',
  'Undifferentiated sarcoma', 'US', 'Undifferentiated sarcomas',
  'Control (reactive tissue)', 'CTRL1', 'Non-neoplastic',
  'Control (blood)', 'CTRL2', 'Non-neoplastic',
  'Control (muscle tissue)', 'CTRL3', 'Non-neoplastic',
  'Chordoma (de-differentiated)', 'CHORD-DD', 'Notochordal tumors',
  'Lipoma', 'LPO', 'Adipocytic tumors',
  'Atypical lipomatous tumour', 'ALT/WDLPS', 'Adipocytic tumors',
  'Dermatofibrosarcoma protuberans', 'DFSP', 'Fibroblastic/myofibroblastic tumors',
  'Giant cell tumour of bone', 'GCTB', 'Osteoclastic giant cell rich tumors',
  'Sarcoma not otherwise specified', 'SARC-NOS', 'Undifferentiated sarcomas',
  'Osteoblastoma', 'OBL', 'Osteogenic tumors',
  'Desmoid-type fibromatosis', 'DES', 'Fibroblastic/myofibroblastic tumors',
  'Low-grade fibromyxoid sarcoma', 'LGFMS', 'Fibroblastic/myofibroblastic tumors',
  'Malignant tumour not otherwise specified', 'MAL-NOS', 'Undifferentiated sarcomas',
  'Neurofibroma', 'NF', 'Peripheral nerve sheath tumors',
  'Neurofibroma (plexiform)', 'NF-P', 'Peripheral nerve sheath tumors',
  'Schwannoma', 'SCHW', 'Peripheral nerve sheath tumors',
  'Infantile fibrosarcoma', 'IFS', 'Fibroblastic/myofibroblastic tumors',
  'Kaposi sarcoma', 'KS', 'Vascular tumors',
  'Melanoma', 'MEL', 'Melanocytic tumors',
  'Squamous cell carcinoma (cutaneous)', 'SCC', 'Epithelial tumors',
  'Sclerosing epithelioid sarcoma', 'SEF', 'Fibroblastic/myofibroblastic tumors',
  'Myxofibrosarcoma', 'MFS', 'Fibroblastic/myofibroblastic tumors',
  'Pleomorphic liposarcoma', 'PLPS', 'Adipocytic tumors',
  'Dedifferentiated liposarcoma', 'DDLPS', 'Adipocytic tumors',
  'Well-differentiated liposarcoma', 'ALT/WDLPS', 'Adipocytic tumors',
  'Endometrial stromal sarcoma (high-grade)', 'ESS-HG', 'Endometrial stromal tumors',
  'Endometrial stromal sarcoma (low-grade)', 'ESS-LG', 'Endometrial stromal tumors',
  'Epithelioid haemangioendothelioma', 'EHE', 'Vascular tumors',
  'Epithelioid sarcoma', 'ES', 'Tumors of uncertain differentiation',
  'Epithelioid sarcoma, classic', 'ES', 'Tumors of uncertain differentiation',
  'Langerhans cell histiocytosis', 'LCH', 'Histiocytic neoplasms',
  'Myositis ossificans', 'MOS', 'Fibroblastic/myofibroblastic tumors',
  'Myositis proliferans', 'MYP', 'Fibroblastic/myofibroblastic tumors',
  'Nodular fasciitis', 'NOD', 'Fibroblastic/myofibroblastic tumors',
  'Ossifying fibromyxoid tumour', 'OFMT', 'Tumors of uncertain differentiation',
  'Glioblastoma', 'GBM', 'Glial tumors',
  'Inflammatory myofibroblastic tumour', 'IMT', 'Fibroblastic/myofibroblastic tumors'
)

# Copy the color pallette from Koelsche et al 2021
# The colors are assigned to the 65 methylation classes, not to unique diagnoses
# Add colors so that each of the diagnoses is assigned a unique color
color_pal <- unique(temp_meta$Color)
unique_dx <- unique(dx_map$Dx)
new_colors <- pals::glasbey(n = length(unique_dx) - length(color_pal))
color_pal <- unique(c(color_pal, new_colors))
message(sprintf("Unique diagnoses: %d", length(unique_dx)))
message(sprintf("Unique colors: %d", length(color_pal)))
names(color_pal) <- unique_dx

dx_map <- dx_map %>%
  mutate(Color_dx = color_pal[Dx])

message("   Generating WHO Colors...")
unique_who <- sort(unique(dx_map$WHO_differentiation))
who_colors <- pals::alphabet(n = length(unique_who))
message(sprintf("Unique WHO categories: %d", length(unique_who)))
message(sprintf("Unique WHO colors: %d", length(who_colors)))
names(who_colors) <- unique_who
dx_map <- dx_map %>% mutate(Color_WHO = who_colors[WHO_differentiation])

# Merge with temp_meta
meta_complete <- temp_meta %>%
    left_join(dx_map, by = "Diagnosis") 
    
# Save
message("Saving metadata file")
write.csv(meta_complete, file.path(processed_dir, "metadata_complete.csv"), row.names = FALSE)

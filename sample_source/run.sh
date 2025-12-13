#!/bin/bash
# filepath: /Users/tahsin/Desktop/ceng477-hw2/sample_source/run.sh

# Renk kodları
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Ana çıktı klasörü
OUTPUT_BASE="my_outputs"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   Rasterizer Build & Test Script${NC}"
echo -e "${BLUE}========================================${NC}"

# Build
echo -e "\n${YELLOW}[BUILD] Derleniyor...${NC}"
make clean
make

if [ $? -ne 0 ]; then
    echo -e "${RED}[HATA] Derleme başarısız!${NC}"
    exit 1
fi

echo -e "${GREEN}[BUILD] Derleme başarılı!${NC}\n"

# Fonksiyon: Input klasörünü çalıştır ve çıktıları taşı
run_and_move() {
    local input_dir="$1"
    local output_subdir="$2"
    local output_dir="$OUTPUT_BASE/$output_subdir"
    
    if [ -d "$input_dir" ]; then
        echo -e "${BLUE}----------------------------------------${NC}"
        echo -e "${BLUE}Input:  $input_dir${NC}"
        echo -e "${BLUE}Output: $output_dir${NC}"
        echo -e "${BLUE}----------------------------------------${NC}"
        
        # Output klasörünü oluştur (yoksa)
        mkdir -p "$output_dir"
        
        for xml_file in "$input_dir"/*.xml; do
            if [ -f "$xml_file" ]; then
                filename=$(basename "$xml_file")
                echo -e "${YELLOW}[ÇALIŞTIRILIYOR] $filename${NC}"
                ./rasterizer "$xml_file"
                
                if [ $? -eq 0 ]; then
                    echo -e "${GREEN}[TAMAM] $filename${NC}"
                else
                    echo -e "${RED}[HATA] $filename${NC}"
                fi
            fi
        done
        
        # Bu input için oluşan PPM dosyalarını ilgili output klasörüne taşı
        if ls *.ppm 1> /dev/null 2>&1; then
            ppm_count=$(ls *.ppm | wc -l)
            mv *.ppm "$output_dir/" 2>/dev/null
            echo -e "${GREEN}[TAŞINDI] $ppm_count PPM dosyası -> $output_dir/${NC}"
        fi
        
        echo ""
    fi
}

# Her input/output çiftini çalıştır
run_and_move "inputs_outputs/culling_disabled_inputs" "culling_disabled_outputs"
run_and_move "inputs_outputs/culling_enabled_inputs" "culling_enabled_outputs"
run_and_move "inputs_outputs/clipping_example" "clipping_example"
run_and_move "inputs_outputs/different_projection_type" "different_projection_type"

echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}   Tüm testler tamamlandı!${NC}"
echo -e "${BLUE}========================================${NC}"

# Tüm output klasörlerindeki dosyaları listele
echo -e "\n${BLUE}Oluşturulan çıktılar:${NC}"
for output_subdir in \
    "culling_disabled_outputs" \
    "culling_enabled_outputs" \
    "clipping_example" \
    "different_projection_type"
do
    output_dir="$OUTPUT_BASE/$output_subdir"
    if [ -d "$output_dir" ]; then
        if ls "$output_dir"/*.ppm 1> /dev/null 2>&1; then
            count=$(ls "$output_dir"/*.ppm | wc -l)
            echo -e "${YELLOW}$output_dir: $count dosya${NC}"
        fi
    fi
done
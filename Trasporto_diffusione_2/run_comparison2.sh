#!/bin/bash

# Script per eseguire le soluzioni stabilizzate e non stabilizzate per diversi valori di 'a'
# Output: file dati per MATLAB

# Colori per output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Valori del parametro 'a' da testare
A_VALUES=(0.130 0.129 0.128 0.127 0.126 0.125 )

# Parametri per le simulazioni
FE_DEGREE=1
INITIAL_REFINEMENT=2
N_CYCLES=6

# Eseguibili
EXEC_NOSTAB="./TrasDiff2"
EXEC_STAB="./TrasDiff3"

# Directory per i risultati
RESULTS_DIR="results_td2vstd3"
mkdir -p $RESULTS_DIR

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Script stabilizzato + non stabilizzato${NC}"
echo -e "${GREEN}========================================${NC}"

# Funzione per creare il file dei parametri
create_param_file() {
    local a_val=$1
    local filename=$2

    cat > $filename << EOF
subsection Problem parameters
  set Finite element degree = $FE_DEGREE
  set Initial refinement = $INITIAL_REFINEMENT
  set Number of cycles = $N_CYCLES
  set Exact solution expression = exp(-(x+y)/$a_val)
  set Right hand side expression = -3*exp(-(x+y)/$a_val)/$a_val
  set Coefficiente diffusione = $a_val
end

subsection Convergence table
  set Enable computation of the errors = true
  set Error file name                  = errors2.txt
  set Error precision                  = 3
  set Exponent for p-norms             = 2
  set Extra columns                    = dofs, cells
  set List of error norms to compute   = Linfty_norm, L2_norm, H1_norm
  set Rate key                         = dofs
  set Rate mode                        = reduction_rate_log2
end

EOF
}

# File di output MATLAB
OUTPUT_FILE="$RESULTS_DIR/tabella_td2_td3.txt"

echo "% Dati per FEM stabilizzato e non stabilizzato" > $OUTPUT_FILE
echo "% Colonne: a_value, peclet, L2_nostab, H1_nostab, gmres_nostab, L2_stab, H1_stab, gmres_stab" >> $OUTPUT_FILE
echo "%" >> $OUTPUT_FILE

# Loop su tutti i valori di 'a'
for a in "${A_VALUES[@]}"; do
    peclet=$(echo "scale=6; 1.0 / (2 * $a)" | bc)

    echo -e "\n${YELLOW}Testando a = $a (Peclet ≈ $peclet)${NC}"

    # Crea file dei parametri
    create_param_file $a "parametri_exp_2d.prm"

    ##############################
    #   NON STABILIZZATO
    ##############################
    echo "  → Eseguendo non stabilizzato..."

    if [ -f "$EXEC_NOSTAB" ]; then
        $EXEC_NOSTAB > $RESULTS_DIR/output_nostab_a${a}.txt 2>&1

        error_l2_nos=$(grep -A 100 "cells" $RESULTS_DIR/output_nostab_a${a}.txt | grep -E "^\s*[0-9]+" | tail -1 | awk '{print $3}')
        error_h1_nos=$(grep -A 100 "cells" $RESULTS_DIR/output_nostab_a${a}.txt | grep -E "^\s*[0-9]+" | tail -1 | awk '{print $5}')
        gmres_nos=$(grep "GMRES iterations needed" $RESULTS_DIR/output_nostab_a${a}.txt | tail -1 | awk '{print $1}')

        echo "    L2: $error_l2_nos, H1: $error_h1_nos, GMRES: $gmres_nos"
    else
        echo "    ERRORE: $EXEC_NOSTAB non trovato!"
        error_l2_nos="NaN"
        error_h1_nos="NaN"
        gmres_nos="NaN"
    fi

    ##############################
    #   STABILIZZATO
    ##############################
    echo "  → Eseguendo stabilizzato..."

    if [ -f "$EXEC_STAB" ]; then
        $EXEC_STAB > $RESULTS_DIR/output_stab_a${a}.txt 2>&1

        error_l2_stab=$(grep -A 100 "cells" $RESULTS_DIR/output_stab_a${a}.txt | grep -E "^\s*[0-9]+" | tail -1 | awk '{print $3}')
        error_h1_stab=$(grep -A 100 "cells" $RESULTS_DIR/output_stab_a${a}.txt | grep -E "^\s*[0-9]+" | tail -1 | awk '{print $5}')
        gmres_stab=$(grep "GMRES iterations needed" $RESULTS_DIR/output_stab_a${a}.txt | tail -1 | awk '{print $1}')

        echo "    L2: $error_l2_stab, H1: $error_h1_stab, GMRES: $gmres_stab"
    else
        echo "    ERRORE: $EXEC_STAB non trovato!"
        error_l2_stab="NaN"
        error_h1_stab="NaN"
        gmres_stab="NaN"
    fi

    # Salvataggio nel file MATLAB
    echo "$a $peclet $error_l2_nos $error_h1_nos $gmres_nos $error_l2_stab $error_h1_stab $gmres_stab" >> $OUTPUT_FILE

done

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}Completato!${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "File dati salvato in: ${YELLOW}$OUTPUT_FILE${NC}"

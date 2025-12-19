#!/bin/bash
# run_convergence_folder.sh
# Script per testare l'ordine di convergenza SUPG al variare di b

# Lista dei valori di b da testare
b_values=(1 10 100 1000 10000 100000 1000000 10000000 100000000 1000000000) 

# Numero di cicli di raffinamento
n_cycles=8

# Degree FE
fe_degree=1

# Eseguibile SUPG
exe_supg="./TrasDiffComparison_supg"

# Cartella dei risultati
results_dir="results_gls_supg_0.025_80-85"
mkdir -p $results_dir

# File di output finale
output_file="$results_dir/convergence_comparison.txt"
echo -e "b\tMethod\tCells\tDOFs\tL2_error\tH1_error\tL2_order\tGMRES_iter" > $output_file

for b in "${b_values[@]}"; do
    echo "==============================="
    echo "Testing b = $b"

    # Creiamo il file di parametri
    param_file="$results_dir/param_b_${b}.prm"
    cat > $param_file <<EOL
    
subsection Problem parameters
  set Finite element degree = $fe_degree
  set Initial refinement = 2
  set Number of cycles = $n_cycles
  set Exact solution expression = exp(-(x+y)/a)
  set Right hand side expression = -(2+2*b)*exp(-(x+y)/a)/a+exp(-(x+y)/a)  
  set Coefficiente diffusione = 0.025
  set Coefficiente trasporto = $b
end
subsection Convergence table
  set Enable computation of the errors = true
  set Error file name                  = $results_dir/errors_b_${b}.txt
  set Error precision                  = 3
  set Exponent for p-norms             = 2
  set Extra columns                    = dofs, cells
  set List of error norms to compute   = L2_norm, H1_norm
  set Rate key                         = dofs
  set Rate mode                        = reduction_rate_log2
end
EOL

    # Copiamo il file nella posizione richiesta dal solver
    cp $param_file parametri_comp_2d.prm

    # --- ESECUZIONE SUPG ---
    echo "Running SUPG..."
    log_file="$results_dir/log_SUPG_b_${b}.txt"
    $exe_supg > $log_file

    # Estrazione errori e ordini
    last_line=$(tail -n 1 "$results_dir/errors_b_${b}.txt")
    cells=$(echo $last_line | awk '{print $1}')
    dofs=$(echo $last_line | awk '{print $2}')
    L2_error=$(echo $last_line | awk '{print $3}')
    L2_order=$(echo $last_line | awk '{print $4}')
    H1_error=$(echo $last_line | awk '{print $5}')

    # Estrazione iterazioni GMRES
    gmres_iter=$(grep "GMRES iterations needed to obtain convergence" $log_file | tail -n 1 | awk '{print $1}')

    # Salvataggio
    echo -e "${b}\tSUPG\t${cells}\t${dofs}\t${L2_error}\t${H1_error}\t${L2_order}\t${gmres_iter}" >> $output_file

done

echo "Confronto completato. Tutti i risultati sono in $results_dir"

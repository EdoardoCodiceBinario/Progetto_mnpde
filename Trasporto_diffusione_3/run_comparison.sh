#!/bin/bash

EXEC="./TrasDiff3"
PARAM_FILE="parametri_exp_2d.prm"

A_VALUES=(1 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15)

OUTDIR="results_exp"
LOGDIR="$OUTDIR/logs"
SUMMARY="$OUTDIR/summary_table.txt"

mkdir -p "$LOGDIR"
> "$SUMMARY"

echo "# a  Peclet  L2_error  H1_error  GMRES_iter  L2_order" >> "$SUMMARY"

create_param_file() {
    local a=$1
    cat > "$PARAM_FILE" << EOF
subsection Problem parameters
  set Finite element degree = 1
  set Initial refinement = 2
  set Number of cycles = 6
  set Exact solution expression = exp(-(x+y)/$a)
  set Right hand side expression = -3*exp(-(x+y)/$a)/$a
  set Coefficiente diffusione = $a
end

subsection Convergence table
  set Enable computation of the errors = true
  set Error file name                  = errors2.txt
  set Error precision                  = 8
  set Exponent for p-norms             = 2
  set List of error norms to compute   = L2_norm, H1_norm
  set Rate key                         = cells
  set Rate mode                        = reduction_rate_log2
end
EOF
}

for a in "${A_VALUES[@]}"; do
    echo "Running a = $a"

    create_param_file "$a"

    LOGFILE="$LOGDIR/log_a_${a}.txt"
    $EXEC > "$LOGFILE" 2>&1

    # ultima riga numerica della tabella
    last_row=$(grep -E "^[[:space:]]*[0-9]+" "$LOGFILE" | tail -n 1)

    L2=$(echo "$last_row" | awk '{print $3}')
    L2_order=$(echo "$last_row" | awk '{print $4}')
    H1=$(echo "$last_row" | awk '{print $5}')

    GMRES=$(grep "GMRES iterations" "$LOGFILE" | tail -n 1 | awk '{print $1}')

    # Peclet (b = 1, h = 2^-6)
    h=$(awk 'BEGIN{print 2^(-6)}')
    Peclet=$(awk -v h="$h" -v a="$a" 'BEGIN{printf "%.3e", h/(2*a)}')

    echo "$a $Peclet $L2 $H1 $GMRES $L2_order" >> "$SUMMARY"
done

echo "Fatto."
echo "Tabella finale: $SUMMARY"
echo "Log completi in: $LOGDIR"

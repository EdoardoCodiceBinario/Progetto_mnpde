% Legge il file di errore e fa il plot log-log delle norme L2 e Linfty
data = readmatrix('errors1_brutto.txt','FileType','text');

% Numero di iterazioni
k = 7;
v = 1:k;

% Estrai le colonne
L2_norm = data(:,3);
Linfty_norm = data(:,5);

% Rimuovi eventuali NaN
mask_L2 = ~isnan(L2_norm);
mask_Linf = ~isnan(Linfty_norm);

% Plot log-log
figure;
hold on
plot(v(mask_L2), L2_norm(mask_L2), '-o', 'LineWidth', 1.5);
plot(v(mask_Linf), Linfty_norm(mask_Linf), '-s', 'LineWidth', 1.5);

xlabel('Iterazione / Numero di celle');
ylabel('Errore');
title('Errore di approssimazione');
grid on
legend('||u - u_h||_{L2}','||u - u_h||_{\infty}','Location','northeast');
hold off

% === Salvataggio del grafico ===
% (il comando saveas salva l'ultima figura aperta)
saveas(gcf, 'plot_td1_linf_l2_brutto.jpg');

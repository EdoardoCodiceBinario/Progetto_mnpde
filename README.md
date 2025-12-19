# Progetto_mnpde
Questa repository contiene i file utilizzati nel corso della realizzazione del progetto per il corso "metodi numerici per le equazioni alle derivate parziali", precisamente riguardo a metodi di stabilizzazione per problemi di trasporto, diffusione e reazione.

La formulazione che è stata considerata è la seguente

$$
\begin{cases}
-\nabla \cdot (a \nabla u) + \mathbf{b} \cdot \nabla u + c u = f & \text{in } \Omega, \\
u = 0 & \text{on } \partial\Omega.
\end{cases}
$$

Segue il contenuto delle cartelle 


Cartella TrasDiff2:

Consideriamo il problema 

$$
u(x,y) = e^{-\frac{x+y}{a}}.
$$

L'aperto d'interesse è $\Omega = (-1,1)^2$.  
La formulazione forte del problema è

$$
\begin{cases}
-a \Delta u + \dfrac{\partial u}{\partial x}
= -\dfrac{3}{a} e^{-\frac{x+y}{a}} & \text{in } \Omega, \\
u = u(x,y) & \text{on } \partial\Omega.
\end{cases}
$$

Nella cartella si trova il file e l'eseguibile che risolve il problema senza stabilizzazione, per valori minori di 0.125 
la dominazione per diffusione è troppo pronunciata è quindi non si arriva a convergenza iterativa.
Inoltre, si po


Cartella TrasDiff3:

Il problema è il medesimo del precedente, però l'aperto considerato è $\Omega = (0,1)^2$,
questo per limitare l'importanza dello strato limite nel fenomeno della stabilizzazione, 
per chi fosse interessato a possibili risoluzioni del problema faccio riferimento a [metti riferimento].   


Nella carte


Cartella confronto:




Cartella problema_senza_soluzione_esplicita:


# Progetto_mnpde

Questa repository contiene i file utilizzati durante la realizzazione del progetto per il corso  
**Metodi numerici per le equazioni alle derivate parziali**, con particolare attenzione ai
metodi di stabilizzazione per problemi di **trasporto, diffusione e reazione**.

La formulazione considerata è la seguente:

$$
\begin{cases}
-\nabla \cdot (a \nabla u) + \mathbf{b} \cdot \nabla u + c u = f & \text{in } \Omega, \\
u = 0 & \text{on } \partial\Omega.
\end{cases}
$$

Se si prova a risolvere rispetto ai metodi Galerkin standard il problema si nota che per valori di a che tendono a 0, e per valori
di b che crescono vengono introdotte oscillazioni numeriche nel sistema che portano a non ottenre una convergenza iterativa.

I metodi che vengono dunque proposti sono metodi di Petrov Galerkin detti di stabilizzazione, nel senso di eliminazione delle suddette oscillazioni.

I metodi che sono presentati sono il Galerkin Least Square(GLS) e lo Streamline Upwind Petrov Galerkin(SUPG).

Si tratta di cambiare il problema variazionale in un'altro che tenda al primo al raffinarsi della approssimazione.

Di seguito è riportato il contenuto delle principali cartelle del repository.

---

## **Cartella TrasDiff2**

Consideriamo il problema con soluzione nota

$$
u(x,y) = e^{-\frac{x+y}{a}}.
$$

Dove a è un parametro che andremo a variare per studiare il fenomeno della dominazione per diffusione.
L’aperto di interesse è $\Omega = (-1,1)^2$.  
La formulazione forte del problema è

$$
\begin{cases}
-a \Delta u + \dfrac{\partial u}{\partial x}
= -\dfrac{3}{a} e^{-\frac{x+y}{a}} & \text{in } \Omega, \\
u = u(x,y) & \text{on } \partial\Omega.
\end{cases}
$$

Nella cartella sono presenti il file sorgente e l’eseguibile che risolvono il problema
**senza stabilizzazione**.  
Per valori di $a < 0.125$ la dominazione del termine di trasporto è troppo pronunciata e
non si osserva convergenza del metodo iterativo.

Nella cartella si trova:

TrasDiff2.cpp/TrasDiff2 ----> file con metodo fem classico, metodo Galerkin
parametri_exp_2d.prm -------> file di parametri per il metodo
run_comparison2.sh ---------> file shell per confrontare il metodo stabilizzato, i cui file si trovano nella cartella TrasDiff3 con quello non stabilizzato
results_td2vstd3 -----------> cartella con i risultati dell'errore rispetto alla soluzione originali per valori vicino al punto critico di convergenza iterativa, sono presenti anche gli output da terminale dei due file, l'ultima iterazione è salvata, al variare di a, nel file tabella_td2_td3.txt
---

## **Cartella TrasDiff3**

Il problema considerato è il medesimo del caso precedente, ma l’aperto è

$$
\Omega = (0,1)^2,
$$

al fine di ridurre l’influenza dello **strato limite** nel fenomeno della stabilizzazione.  

Per chi fosse interessato a possibili strategie di risoluzione del problema, si rimanda a
[**riferimento da inserire**].

Mediante stabilizzazione col metodo Galerkin Least Square(GLS) si arriva a convergenza iterativa
per valori arbitrati di a, si osserva che l'ordine di convergenza dell'errore in norma L2 si stabilizza
al valore di 1/2.

---

## **Cartella Confronto**

In questa cartella viene effettuato un confronto tra i metodi di stabilizzazione **SUPG**
e **GLS**, analizzando il comportamento in termini di oscillazioni numeriche ed errori al
variare dei parametri del problema.

---

## **Cartella Problema_senza_soluzione_esplicita**

In questa cartella viene studiato un problema di trasporto-diffusione-reazione per il quale
non è disponibile una soluzione analitica esplicita.  
L’analisi è condotta esclusivamente tramite risultati numerici.

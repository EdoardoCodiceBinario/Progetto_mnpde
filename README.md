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

## **Cartella Trasporto_Diffusione_2**

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


Per valori di $a < 0.125$ la dominazione del termine di trasporto è troppo pronunciata e
non si osserva convergenza iterativa.

Questo esempio evidenzia le problematiche che insorgono quando si va a risolvere mediante metodi FEM Galerkin
problemi a diffusione dominante.



| File / Cartella | Descrizione |
|-----------------|-------------|
| **TrasDiff2.cpp / TrasDiff2** | Implementazione del metodo FEM classico (Galerkin). |
| **parametri_exp_2d.prm** | File di parametri del metodo. |
| **run_comparison2.sh** | Script per il confronto tra metodo stabilizzato (TrasDiff3) e non stabilizzato. |
| **results_td2vstd3/** | Cartella degli output terminale e il file **tabella_td2_td3.txt** riassuntivo dell'ultima iterazione al variare di a. |

---

## **Cartella Trasporto_Diffusione_3**

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

Osserviamo che per la formulazione del problema i metodi SUPG e GLS coincidono.

| File / Cartella | Descrizione |
|-----------------|-------------|
| **TrasDiff3.cpp / TrasDiff3** | Implementazione del metodo GLS/SUPG. |
| **parametri_exp_2d.prm** | File di parametri del metodo. |
| **run_comparison.sh** | Script per il calcolo dell'errore e ordine di convergenza  |
| **results_exp** | Risultati dell’errore per valori di $a$ da 1 a e-15, un riassunto dei valori all'ultima iterazione di raffinamento si trova a **summary_table.txt**. |

Da l'ultima tabella sono stati ricavati inoltre l'ordine di convergenza dell'errore di approssimazione in norma L2, il quale è importante sottlineare che abbiamo ordine di convergenza 1/2.

---

## **Cartella Confronto**

In questa cartella viene effettuato un confronto tra i metodi di stabilizzazione **SUPG**
e **GLS**, analizzando il comportamento in termini di oscillazioni numeriche ed errori al
variare dei parametri del problema.
Il problema studiato ha la medesima soluzione dei precedenti ma una formulazione differente

$$\begin{cases}
-a \Delta u + \mathbf{b}\cdot\nabla u +  u = f & \text{in } \Omega, 
u = u(x,y) & \text{su } \partial \Omega.
\end{cases}$$

dove

$$\mathbf{b}(x,y) = b \cdot
\begin{bmatrix}
1 \\
1
\end{bmatrix},$$

E' stato dunque introdotto un parametro per il trasporto mentre quello di diffusione è stato posto a 0.025, si nota che per valori
di $b>80$ il metodo GLS non converge, mentre il metodo SUPG mostra convergenza più che ottimale, facendo riferimento alla convergenza prevista da Lemma di Ceà+Lemma di Bramble-Hilbert.

| File / Cartella | Descrizione |
|-----------------|-------------|
| **TrasDiffComparison_gls.cpp / TrasDiffComparison_gls** | Implementazione del metodo GLS. |
| **TrasDiffComparison_supg.cpp / TrasDiffComparison_supg** | Implementazione del metodo SUPG. |
| **parametri_comp_2d.prm** | File di parametri dei metodi. |
| **comparison_gls_supg_0.025_1-e+9.sh** | Script per il calcolo dell'errore e ordine di convergenza nel caso del metodo SUPG, i valori considerati sono da b che varia da 1 a e+9 |
| **results_gls_supg_0.025** | Cartella dove si calcola i soliti parametri e l'ordine di convergenza dei due metodi a confronti per valori di b da 1 a 80  |

---

## **Cartella Problema_senza_soluzione_esplicita**

In questa cartella troviamo un esempio di problema che non sappiamo risolvere esplicitamente e nel quale la stabilizzazione è cruciale.
Il problema studiato ha formulazione differenziale

$$
\begin{cases}
-\nabla \cdot (a \nabla u) + \mathbf{b} \cdot \nabla u + c u = f & \text{in } \Omega, \\
u = 0 & \text{su } \partial \Omega.
\end{cases}
$$

dove abbiamo:

$$
\mathbf{b}(x,y) =
\begin{bmatrix}
30 (y-0.5) + 50 \sin(5 \pi y) \\
-30 (x-0.5) + 50 \cos(5 \pi x)
\end{bmatrix},
\quad
c(x,y) = 50 \, e^{-50((x-0.5)^2 + (y-0.5)^2)} + 15 \, \sin(5 \pi x) \cos(5 \pi y) \\ 
f(x,y) = 20 \, \exp\Big(- (x-0.3)^2 - (y-0.1)^2 \Big)
$$


| File / Cartella | Descrizione |
|-----------------|-------------|
| **TrasDiffboh.cpp / TrasDiffboh** | Implementazione del metodo classico. |
| **TrasDiffboh_gls.cpp / TrasDiffboh_gls** | Implementazione del metodo GLS. |
| **parametri_boh_2d.prm** | File di parametri. |
---

## **Cartella Oscillazioni**

In questa cartella si produce un esempio dove si vede esplicitamente il fenomeno oscillatorio, combinato con uno strato limite presente
sul bordo nelle vicinanze del punto (1,1).

Andiamo a risolvere un altro problema:

$$
\begin{cases}
-a \Delta u + \dfrac{\partial u}{\partial x} + \dfrac{\partial u}{\partial y} = 1 & \text{in } \Omega, \\
u(x,y) = 0 & \text{su } \partial\Omega.
\end{cases}
$$

Dove le condizioni sono:


| File / Cartella | Descrizione |
|-----------------|-------------|
| **oscillazioni.cpp / oscillazioni** | Implementazione del metodo classico. |
| **oscillazioni_gls.cpp / oscillazioni_gls** | Implementazione del metodo GLS. |
| **parametri_comp_2d.prm** | File di parametri. |

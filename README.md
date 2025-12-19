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

Di seguito è riportato il contenuto delle principali cartelle del repository.

---

## **Cartella TrasDiff2**

Consideriamo il problema con soluzione nota

$$
u(x,y) = e^{-\frac{x+y}{a}}.
$$

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

---

## **Cartella TrasDiff3**

Il problema considerato è il medesimo del caso precedente, ma l’aperto è

$$
\Omega = (0,1)^2,
$$

al fine di ridurre l’influenza dello **strato limite** nel fenomeno della stabilizzazione.  

Per chi fosse interessato a possibili strategie di risoluzione del problema, si rimanda a
[**riferimento da inserire**].

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

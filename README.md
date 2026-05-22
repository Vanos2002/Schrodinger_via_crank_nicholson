# 1D Časově závislá Schrödingerova rovnice

Numerické řešení 1D časově závislé Schrödingerovy rovnice metodou Crank-Nicolson s vizualizací výsledků v Pythonu.

## Fyzikální popis

Řešíme rovnici

$$i\partial_t \psi(x,t) = \left[-\frac{1}{2}\partial_x^2 + V(x)\right]\psi(x,t)$$

na intervalu $x \in [-25, 25]$, $t \in [0, 5]$, s počáteční podmínkou

$$\psi(x, 0) = \left(\frac{2}{\pi}\right)^{1/4} \exp\left[-(x-1)^2\right]$$

a potenciálem (dvojitá bariéra s absorpčními okraji)

$$V(x) = (x^2 - 1)e^{-x^2/4} - i\cdot 10^{-3}(|x|-8)^4\\Theta(|x|-8)$$

Studuje se únik částice z potenciálové bariéry přes zákon zachování pravděpodobnosti:

$$1 = p(T) + J_+(T) + J_-(T)$$

kde $p(T) = \int_{-a}^{a}|\psi|^2dx$ a $J_\pm(T) = \pm\int_0^T \mathrm{Im}(\psi^*\partial_x\psi)\big|_{x=\pm a}dt$.

## Soubory

| Soubor | Popis |
|---|---|
| `main.cpp` | Simulace – Crank-Nicolson, generování dat |
| `plotting.py` | Vizualizace – vlnové funkce, toky, konvergence |

## Požadavky

**C++**
- kompilátor s podporou C++17 (g++, clang++)
- standardní knihovna (žádné externí závislosti)

**Python**
- Python 3.8+
- `numpy`
- `matplotlib`

```bash
pip install numpy matplotlib
```

## Kompilace a spuštění

```bash
# Kompilace
g++ -O2 -std=c++17 -o main main.cpp

# Spuštění simulace (generuje datové soubory)
./main

# Vizualizace
python3 plotting.py
```

## Výstupní soubory

### Datové soubory (generované C++ programem)

| Vzor názvu | Obsah |
|---|---|
| `wave_t_<T>.dat` | $x$, $\|\psi(x,T)\|^2$ pro referenční síť |
| `free_t_<T>.dat` | $x$, $\|\psi(x,T)\|^2$ pro volnou částici ($V=0$) |
| `wave_prob.dat` | $t$, $p$, $j_+$, $j_-$, $J_+$, $J_-$, $p+J_++J_-$ |
| `conv_dx<N>_prob.dat` | Totéž pro různá $\Delta x$ (konvergenční studie) |

### Obrázky (generované Python skriptem)

| Soubor | Obsah |
|---|---|
| `fig0_potential.png` | Tvar potenciálu a počáteční vlnová funkce |
| `fig1_wave_functions.png` | $\|\psi(x,t)\|^2$ pro $t = 0, 1, 2, 3, 4, 5$ + srovnání s volnou částicí |
| `fig2_time_evolution.png` | Časový vývoj $p(T)$, $J_\pm(T)$, odchylka od zachování |
| `fig3_convergence.png` | Konvergence $p$, $J_\pm$ při $T=1$, $a=7$ pro různá $\Delta x$ |

## Numerická metoda

### Diskretizace

| Parametr | Referenční hodnota |
|---|---|
| $\Delta x$ | $0{,}05$ |
| $\Delta t$ | $0{,}002$ |
| $N$ (počet bodů) | $1001$ |

Crank-Nicolsonovo schéma je **bezpodmínečně stabilní** a má přesnost $O(\Delta x^2, \Delta t^2)$. Výsledná tridiagonální soustava se řeší Thomasovým algoritmem v čase $O(N)$.

Okrajové podmínky jsou Dirichletovy ($\psi = 0$ na krajích), absorpční imaginární část potenciálu pro $|x| > 8$ zabraňuje odrazu od hranic domény.

### Konvergenční studie

Program automaticky spustí tři simulace s $\Delta x \in \{0{,}1,\; 0{,}05,\; 0{,}025\}$ a $\Delta t = \Delta x^2 / 4$ do času $T = 1$. Výsledky se porovnávají s nejjemnější sítí.

## Výsledky

Při $T = 1$, $a = 7$ je prakticky celá pravděpodobnost stále uvnitř oblasti – tunelovací únik začíná být patrný přibližně od $t \approx 1{,}8$ (pravá hranice) a $t \approx 2{,}6$ (levá hranice). Asymetrie plyne z počáteční polohy balíku u $x = 1$, blíže pravé jámě potenciálu.

Zákon zachování $p + J_+ + J_- \approx 1$ je splněn s odchylkou řádu $10^{-3}$, která pochází z absorpce na okrajích domény a numerické chyby druhého řádu.

<img width="1330" height="730" alt="fig0_potential" src="https://github.com/user-attachments/assets/aa7e6d3c-680b-40ac-832a-3264749dc5ee" />

<img width="2080" height="1181" alt="fig1_wave_functions" src="https://github.com/user-attachments/assets/d40d6c2c-234b-4b38-a78f-2aa9dc9dedf4" />

<img width="1488" height="1622" alt="fig2_time_evolution" src="https://github.com/user-attachments/assets/4d2caa0d-8d7f-40ad-96b1-317cd2687154" />

<img width="1774" height="740" alt="fig3_convergence" src="https://github.com/user-attachments/assets/22febcd1-374b-4798-a632-25e906465578" />


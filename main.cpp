#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

// ======================================================
// 1D Časově závislá Schrödingerova rovnice
// Crank-Nicolsonova metoda
// ======================================================
// konstanty a parametry sítě
const double X_MIN = -25.0;
const double X_MAX  =  25.0;
const double PI = 3.14159265358979323846;
const std::complex<double> I(0.0, 1.0);

// ======================================================
// Parametry sítě – předány za běhu pro konvergenční studii
// ======================================================

struct GridParams {
    double dx;
    double dt;
    int    N;
};

GridParams makeGrid(double dx, double dt) {
    GridParams g;
    g.dx = dx;
    g.dt = dt;
    g.N  = static_cast<int>((X_MAX - X_MIN) / dx) + 1;
    return g;
}

// ======================================================
// Počáteční vlnová funkce  ψ(x,0) = (2/π)^(1/4) exp[-(x-1)²]
// ======================================================

std::vector<std::complex<double>> initializeWaveFunction(const GridParams& g) {
    std::vector<std::complex<double>> psi(g.N);
    const double norm = std::pow(2.0 / PI, 0.25);
    for (int i = 0; i < g.N; ++i) {
        double x = X_MIN + i * g.dx;
        psi[i] = norm * std::exp(-(x - 1.0) * (x - 1.0));
    }
    return psi;
}

// ======================================================
// Volná částice (V=0) – pro srovnání
// ======================================================

std::vector<std::complex<double>> initializeFreePotential(const GridParams& g) {
    return std::vector<std::complex<double>>(g.N, 0.0);
}

// ======================================================
// Potenciál: reálná část "double-well" + imaginární absorpce pro |x|>8
// ======================================================

std::vector<std::complex<double>> initializePotential(const GridParams& g) {
    std::vector<std::complex<double>> V(g.N);
    for (int i = 0; i < g.N; ++i) {
        double x = X_MIN + i * g.dx;
        std::complex<double> V_real = (x * x - 1.0) * std::exp(-x * x / 4.0);
        std::complex<double> V_abs  = 0.0;
        if (std::abs(x) > 8.0) {
            double d = std::abs(x) - 8.0;
            V_abs = -I * 1e-3 * d * d * d * d;
        }
        V[i] = V_real + V_abs;
    }
    return V;
}

// ======================================================
// Crank-Nicolsonova metoda pro jeden časový krok: A ψ_new = B ψ_old
// ======================================================

void crankNicolsonStep(
    std::vector<std::complex<double>>& psi,
    const std::vector<std::complex<double>>& V,
    const GridParams& g)
{
    int N = g.N;
    double dx = g.dx, dt = g.dt;

    std::vector<std::complex<double>> a(N), b(N), c(N), d(N);
    std::vector<std::complex<double>> alpha(N), beta(N);
    std::vector<std::complex<double>> psi_new(N, 0.0);

    std::complex<double> r = I * dt / (2.0 * dx * dx);   // = iΔt/(2Δx²)

    // ---------- Koeficienty levé strany (matice A) ----------
    for (int i = 0; i < N; ++i) {
        a[i] = -r / 2.0;                              // sub-diagonála (pod hlavní diagonálou)
        b[i] =  1.0 + r + I * dt * V[i] / 2.0;       // diagonála
        c[i] = -r / 2.0;                              // super-diagonála (nad hlavní diagonálou)
    }

    // Dirichletovy okrajové podmínky ψ=0
    b[0]   = 1.0; c[0]   = 0.0;
    b[N-1] = 1.0; a[N-1] = 0.0;
    d[0]   = 0.0; d[N-1] = 0.0;

    // ---------- Pravá strana (matice B · ψ) ----------
    for (int i = 1; i < N - 1; ++i) {
        d[i] =   (r / 2.0) * psi[i-1]
               + (1.0 - r - I * dt * V[i] / 2.0) * psi[i]
               + (r / 2.0) * psi[i+1];
    }

    // ---------- Dopředná eliminace ----------
    alpha[0] = c[0] / b[0];
    beta[0]  = d[0] / b[0];
    for (int i = 1; i < N; ++i) {
        std::complex<double> denom = b[i] - a[i] * alpha[i-1];
        alpha[i] = c[i] / denom;
        beta[i]  = (d[i] - a[i] * beta[i-1]) / denom;
    }

    // ---------- Zpětná substituce ----------
    psi_new[N-1] = beta[N-1];
    for (int i = N - 2; i >= 0; --i)
        psi_new[i] = beta[i] - alpha[i] * psi_new[i+1];

    psi_new[0]   = 0.0;
    psi_new[N-1] = 0.0;

    psi = psi_new;
}

// ======================================================
// Pravděpodobnost  p = ∫_{-a}^{a} |ψ|² dx
// ======================================================

double calculateProbability(
    const std::vector<std::complex<double>>& psi,
    double a, const GridParams& g)
{
    double p = 0.0;
    for (int i = 0; i < g.N; ++i) {
        double x = X_MIN + i * g.dx;
        if (x >= -a && x <= a)
            p += std::norm(psi[i]) * g.dx;
    }
    return p;
}

// ======================================================
// Okamžitý tok pravděpodobnosti  j = Im(ψ* ∂_x ψ)
// ======================================================

double calculateCurrent(
    const std::vector<std::complex<double>>& psi,
    double x_position, const GridParams& g)
{
    int idx = static_cast<int>((x_position - X_MIN) / g.dx);
    if (idx <= 0 || idx >= g.N - 1) return 0.0;
    std::complex<double> deriv = (psi[idx+1] - psi[idx-1]) / (2.0 * g.dx);
    return std::imag(std::conj(psi[idx]) * deriv);
}

// ======================================================
// Uložení |ψ(x,t)|² do souboru  (pro Python)
// ======================================================

void saveWaveFunction(
    const std::vector<std::complex<double>>& psi,
    double t, const GridParams& g,
    const std::string& prefix)
{
    int fileIndex = static_cast<int>(std::round(1000.0 * t));
    std::string fname = prefix + "_t_" + std::to_string(fileIndex) + ".dat";
    std::ofstream file(fname);
    for (int i = 0; i < g.N; ++i) {
        double x = X_MIN + i * g.dx;
        file << std::setprecision(12) << x << " " << std::norm(psi[i]) << "\n";
    }
}

// ======================================================
// Hlavní simulace pro dané parametry sítě
// ======================================================

void runSimulation(
    double dx, double dt, double T_MAX,
    const std::string& filePrefix,
    bool saveWaves)
{
    GridParams g = makeGrid(dx, dt);

    auto psi  = initializeWaveFunction(g);
    auto V    = initializePotential(g);

    std::ofstream outP(filePrefix + "_prob.dat");
    outP << "t p J_plus J_minus J_plus_cum J_minus_cum conservation\n";

    double J_plus_cumul  = 0.0;   // kumulativní integrál  ∫ j(+a) dt
    double J_minus_cumul = 0.0;   // kumulativní integrál  ∫ j(-a) dt

    // Uložení počátečního stavu
    if (saveWaves) saveWaveFunction(psi, 0.0, g, filePrefix);

    int step = 0;
    for (double t = dt; t <= T_MAX + 1e-12; t += dt, ++step) {

        crankNicolsonStep(psi, V, g);

        double j_plus  =  calculateCurrent(psi,  7.0, g);
        double j_minus = -calculateCurrent(psi, -7.0, g);   // tok ven vlevo = -j(-7)

        J_plus_cumul  += j_plus  * dt;
        J_minus_cumul += j_minus * dt;

        double p = calculateProbability(psi, 7.0, g);

        // zákon zachování: p + J+ + J- ≈ 1
        double conservation = p + J_plus_cumul + J_minus_cumul;

        outP << std::fixed << std::setprecision(8)
             << t           << " "
             << p           << " "
             << j_plus      << " "
             << j_minus     << " "
             << J_plus_cumul  << " "
             << J_minus_cumul << " "
             << conservation << "\n";

        // Ukládáme vlnovou funkci v celých časech 0,1,2,...
        if (saveWaves) {
            double t_round = std::round(t);
            if (std::fabs(t - t_round) < dt / 2.0 && t_round >= 1.0)
                saveWaveFunction(psi, t_round, g, filePrefix);
        }
    }

    std::cout << "Hotovo: " << filePrefix << "  (dx=" << dx << ", dt=" << dt << ")\n";
}

// ======================================================
// Volná částice pro srovnání (analytické řešení je gaussovský balík)
// ======================================================

void runFreeSimulation(double dx, double dt, double T_MAX) {
    GridParams g = makeGrid(dx, dt);
    auto psi = initializeWaveFunction(g);
    auto V   = initializeFreePotential(g);

    for (double t = dt; t <= T_MAX + 1e-12; t += dt) {
        crankNicolsonStep(psi, V, g);
        double t_round = std::round(t);
        if (std::fabs(t - t_round) < dt / 2.0)
            saveWaveFunction(psi, t_round, g, "free");
    }
    std::cout << "Hotovo: volná částice\n";
}

// ======================================================
// main
// ======================================================

int main() {
    const double T_MAX = 5.0;

    // --- Referenční simulace (uložení vln pro vizualizaci) ---
    runSimulation(0.05, 0.002, T_MAX, "wave", true);

    // --- Volná částice pro srovnání ---
    runFreeSimulation(0.05, 0.002, T_MAX);

    // --- Konvergenční studie: různá dx při dt=dx²/4 ---
    std::vector<double> dxs = {0.1, 0.05, 0.025};
    for (double dx : dxs) {
        double dt = dx * dx / 4.0;
        std::string prefix = "conv_dx" + std::to_string((int)(dx * 1000));
        runSimulation(dx, dt, 1.0, prefix, false);   // jen do T=1 pro konvergenční studii
    }

    std::cout << "\nVšechna data vygenerována. Spusťte python_plot.py\n";
    return 0;
}

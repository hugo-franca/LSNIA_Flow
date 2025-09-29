#include "EVPT_Equacoes.h"

/// Essa funcao realiza as iteracoes do metodo de newton pra calcular strain rate
/// Estou assumindo que a funcao regularizadora eh aquela xi_3
double FlowCurve_MetodoDeNewton(MALHA M, double GammaAtual, double Stress)
{
    double stop = 1e+10;

    /// In this case it's a linear equation, pretty easy solution
    if( Stress<M.yield_stress )
        return Stress/M.eta_0;

    /// In the case Stress>=Yield, we use Newton's method for the nonlinear solution
    int iteration = 0;
    while( stop>1e-5 ) {
        double f = M.yield_stress + M.K*pow(GammaAtual, M.power_n) + M.eta_inf*GammaAtual - Stress;
        double f_derivada = M.power_n*M.K*pow(GammaAtual, M.power_n-1.0) + M.eta_inf;
        double GammaNovo = GammaAtual - f/f_derivada;

        PrintDebug("%d: GAMMA ATUAL e Novo: %lf %lf .... %lf %lf .... %lf\n", iteration, GammaAtual, GammaNovo, f, f_derivada, Stress);
        stop = fabs(GammaNovo - GammaAtual);
        GammaAtual = GammaNovo;
        iteration++;
        
        if( iteration>5000 )
            DeuProblema("FlowCurve: Newtons method diverged...\n");
    }
    
    return GammaAtual;
}

/// Essa funcao faz a evolucao do parametro lambda no tempo (structure parameter)
void Evolucao_StructureParameter(MALHA M, double **U, double **V, double **P,
                                double **Txx, double **Txy, double **Tyy,
                                double **EVPT_Lambda_Velho, double **EVPT_Lambda_Novo)
{
    int linha;

    for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
        int i = M.indicesP[linha].i;
        int j = M.indicesP[linha].j;

        // double sigma_dev = EVPT_CalculaSigmaDev(M, i, j, U, V, P, Txx, Txy, Tyy);
        double sigma_dev = 0.0;
        
        // Calculating the equilibrium strain_rate from the flow_curve
        // Note that if sigma_dev<yield_stress, then the strain_rate is always fixed eta_0/sigma_dev
        double strain_rate = FlowCurve_MetodoDeNewton(M, EVPT_Lambda_Velho[i][j], sigma_dev);
        
        // Calculating the equilibrium viscosity. (and doing a little safety check just to avoid division by zero edge case when stress=0)
        // Note that if sigma_dev<yield_stress, this formula will always result in eta_eq = eta_0
        double eta_eq = (fabs(sigma_dev)<1e-10) ? M.eta_0 : sigma_dev/strain_rate;

        // Calculating the equilibrium structure parameter
        // Note that if sigma_dev<yield_stress, this formula will always result in lambda_eq = 1
        double lambda_eq = ( log(eta_eq) - log(M.eta_inf) )/( log(M.eta_0) - log(M.eta_inf) );
        
        if( isnan(sigma_dev) || isnan(strain_rate) || isnan(eta_eq) || isnan(lambda_eq) )
            DeuProblema("PROBLEMA (%d %d): %lf %lf %lf %lf\n", i, j, sigma_dev, strain_rate, eta_eq, lambda_eq);

        // If the material is NOT thixotropic, the structure just immediately goes to the equilibrium state
        if( (M.tipo_modelo!=MODELO_EVPT) || (M.time_eq<1e-5) ) {
            EVPT_Lambda_Novo[i][j] = lambda_eq;

            if( EVPT_Lambda_Novo[i][j]<0.0 /*|| EVPT_Lambda_Novo[i][j]>1.0*/ )
                DeuProblema("aaaaa... Structure fora do alcance (%d, %d): %e %e %e %e\n", i, j, EVPT_Lambda_Novo[i][j], EVPT_Lambda_Velho[i][j], lambda_eq, M.time_eq);
        }
        // If the material is thixotropic, then we have to update the structure with the equation below
        else {
            double convectivo = Conv_Tensor(EVPT_Lambda_Velho, U, V, i, j, "EVPT_Structure", M);
            EVPT_Lambda_Novo[i][j] = EVPT_Lambda_Velho[i][j] + M.dt*( - convectivo + (1.0 - EVPT_Lambda_Velho[i][j]/lambda_eq)/M.time_eq );
        }

        if( EVPT_Lambda_Novo[i][j]>1.0 )
            EVPT_Lambda_Novo[i][j] = 1.0;
        if( EVPT_Lambda_Novo[i][j]<0.0 /*|| EVPT_Lambda_Novo[i][j]>1.0*/ )
            DeuProblema("Structure fora do alcance (%d, %d): %lf %lf %lf %lf\n", i, j, EVPT_Lambda_Novo[i][j], EVPT_Lambda_Velho[i][j], lambda_eq, M.time_eq);
    }



    return;
}
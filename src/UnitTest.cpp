#include<complex>
#include<vector>
#include"HDistribution.h"
#include"VanillaContract.h"
#include"HSimulation.h"
typedef std::numeric_limits< double > dbl;
using namespace std::complex_literals;

std::vector<std::complex<double>> GetCuiGradient(HParams p, std::complex<double> u, double T);
std::vector<std::complex<double>> GetCuiGradient(HParams p, std::complex<double> u, double T)
{
    auto const& [v0, v_bar, rho,kappa, sigma] = p;
    double sigma_times_rho=sigma*rho;
    double sigma_squared=sigma*sigma;
    double half_T=T*.5;
    double two_kappa_v_bar_over_sigma_squared=2*kappa*v_bar/(sigma*sigma);
    double kappa_v_bar_rho_T_over_sigma=kappa*v_bar*rho*T/sigma;
    
    std::complex<double> ui = 1i * u;
    std::complex<double> u_squared = u * u;

    std::complex<double> xi = kappa - sigma_times_rho * ui; // xi = kappa - sigma * rho * u * i
    std::complex<double> xi2 = xi * xi;
    std::complex<double> m = ui + u_squared; // m = u * i + u^2;
    std::complex<double> d = std::sqrt(xi2 + m * sigma_squared); // d = sqrt(pow(xi,2) + m*pow(sigma,2));

    // alp, calp, salp
    std::complex<double> alpha = d * half_T;
    std::complex<double> cosh_alpha = cosh(alpha);
    std::complex<double> sinh_alpha = sinh(alpha);
    std::complex<double> A2_times_v0 = d * cosh_alpha + xi * sinh_alpha;
    std::complex<double> A1 = m * sinh_alpha;
    std::complex<double> A_over_v0 = A1 / A2_times_v0;

    std::complex<double> D = std::log(d) + (kappa - d) * half_T - std::log((d + xi) * 0.5 + (d - xi) * 0.5 * exp(-d * T));
    
//        std::complex<double> g = exp(minus_kappa_v_bar_rho_T_over_sigma * ui);

    // F = S * e^((r - q) * T)
    // characteristic function: y1 = exp(i * log(F / S) * u) * exp(-A + 2 * kappa * b / pow(sigma, 2) * D) * g
    // But we only care about the second exponential, the rest depends only on market parameters and will be computed separately.
    std::complex<double> char_u =std::exp(-v0 * A_over_v0 + two_kappa_v_bar_over_sigma_squared * D - kappa_v_bar_rho_T_over_sigma * ui);

    // B = d * exp(kappa * T / 2) / (A2 * v0);
    double exp_kappa_times_half_T = exp(kappa * half_T); // exp(kappa * T / 2)
    std::complex<double> B = d * exp_kappa_times_half_T / A2_times_v0;
    //std::cout<<"2xi: "<<xi<<"\td: "<<d<<"\tA_1: "<<A1<<"\tA_2: "<<A2_times_v0/v0<<"\tA"<<A_over_v0*v0<<"\tB: "<<B<<"\tD: "<<D<<"\tchf: "<<char_u<<'\n';
    // g = exp(-kappa * b * rho * T * u1 * i / sigma);
    double kappa_v_bar_rho_T = kappa * v_bar * rho * T;  // TAG: PRECOMPUTE
    double minus_kappa_v_bar_rho_T_over_sigma = -kappa_v_bar_rho_T / sigma;  // TAG: PRECOMPUTE

    std::complex<double> H = xi * cosh_alpha + d * sinh_alpha;

    // lnB = log(B);
    std::complex<double> lnB = D;

    // partial b: y3 = y1*(2*kappa*lnB/pow(sigma,2)-kappa*rho*T*u1*i/sigma);
    double two_kappa_over_sigma_squared = two_kappa_v_bar_over_sigma_squared / v_bar;
    double minus_kappa_rho_T_over_sigma = minus_kappa_v_bar_rho_T_over_sigma / v_bar;

    std::complex<double> h_v_bar = two_kappa_over_sigma_squared * lnB + minus_kappa_rho_T_over_sigma * ui;

    // partial rho:
    double minus_kappa_v_bar_t_over_sigma = minus_kappa_v_bar_rho_T_over_sigma / rho; //-kappa * v_bar * T/sigma;

    std::complex<double> sigma_ui_over_d = sigma * ui / d;
    std::complex<double> pd_prho = -xi * sigma_ui_over_d;
    std::complex<double> pA1_prho = m * cosh_alpha * half_T * pd_prho;
    std::complex<double> pA2_prho = -sigma_ui_over_d * H * (1.0 + xi * half_T);
    std::complex<double> pA_prho = (pA1_prho - A_over_v0 * pA2_prho) / A2_times_v0;
    std::complex<double> pd_phrho_minus_pA2_prho_times_d_over_A2 = pd_prho - pA2_prho * d / A2_times_v0;
    std::complex<double> pB_prho = exp_kappa_times_half_T / A2_times_v0 * pd_phrho_minus_pA2_prho_times_d_over_A2;
    std::complex<double> h_rho = -v0 * pA_prho + two_kappa_v_bar_over_sigma_squared * pd_phrho_minus_pA2_prho_times_d_over_A2 / d + minus_kappa_v_bar_t_over_sigma * ui;

    // partial kappa:
    double v_bar_rho_T_over_sigma = v_bar * rho * T / sigma;
    double two_v_bar_over_sigma_squared = two_kappa_v_bar_over_sigma_squared / kappa; // 2 * v_bar / sigma_squared;

    std::complex<double> minus_one_over_sigma_ui = -1.0 / (sigma * ui);
    std::complex<double> pB_pa = minus_one_over_sigma_ui * pB_prho + B * half_T;
    std::complex<double> h_kappa = -v0 * pA_prho * minus_one_over_sigma_ui + two_v_bar_over_sigma_squared * lnB + kappa * two_v_bar_over_sigma_squared * pB_pa / B - v_bar_rho_T_over_sigma * ui;

    // partial sigma:
    double rho_over_sigma = rho / sigma;
    double four_kappa_v_bar_over_sigma_cubed = 4 * kappa * v_bar / pow(sigma, 3);
    double kappa_v_bar_rho_T_over_sigma_squared = kappa_v_bar_rho_T / sigma_squared;
    std::complex<double> pd_pc = (rho_over_sigma - 1.0 / xi) * pd_prho + sigma * u_squared / d;
    std::complex<double> pA1_pc = m * cosh_alpha * half_T * pd_pc;
    std::complex<double> pA2_pc = rho_over_sigma * pA2_prho - 1.0 / ui * (2.0 / (T * xi) + 1.0) * pA1_prho + sigma * half_T * A1;
    std::complex<double> pA_pc = pA1_pc / A2_times_v0 - A_over_v0 / A2_times_v0 * pA2_pc;
    std::complex<double> h_sigma = -v0 * pA_pc - four_kappa_v_bar_over_sigma_cubed * lnB + two_kappa_v_bar_over_sigma_squared / d * (pd_pc - d / A2_times_v0 * pA2_pc) + kappa_v_bar_rho_T_over_sigma_squared * ui;
    return { -A_over_v0 * char_u, char_u * h_v_bar, char_u * h_rho, char_u * h_kappa,  char_u * h_sigma};
}
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"HCalibration.h"
#include"BSM.h"
#include"RNG.h"

#include"levmar.h"
#include<memory.h>
//#define NDEBUG
#define DISTR_TOLERANCE 1e-6
#include<assert.h>
void __assert_fail(const char * assertion, const char * file, unsigned int line, const char * function)
{
    std::cout<<"assertion "<<assertion<<" in file:"<<file<<", function: "<<function<<", line: "<<line<<" failed"<<std::endl;
    abort();
}

void distr_test(){
    double v0 =0.0727878;
    double	v_bar= 0.0946025;
    double rho= 0.0576755;
    double kappa= 1.04079;
    double sigma=1.64704;
    /*
    for(int i=1;i<1000; i++){
        test.p.v_0=v0 * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        test.p.v_m=v_bar * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        test.p.rho=rho * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        test.p.kappa=kappa * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        test.p.sigma=sigma * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        double tau=(((double) rand() / (RAND_MAX)));// *trading_days;
        double x=(((double) rand() / (RAND_MAX)) - 0.5)*100;
        std::complex<double> correct =chf2(test.p,x,tau);
        std::complex<double> wrong = test.chf(x,tau);
        std::cout<<"v_0: "<<test.p.v_0<<"\tv_m: "<<test.p.v_m<<"\trho: "<<test.p.rho<<"\tkappa: "<<test.p.kappa<<"\tsigma: "<<test.p.sigma<<"\ttau: "<<tau<<"\tx: "<<x<<"\tmy_chf: "<<wrong<<"\tEudald Romo chf: "<< correct<<"\tdiff: "<<std::fabs(correct-wrong)<<'\n';
    }*/
    //bool test_passed=true;
    for(int i=1;i<1000; i++){
         HDistribution test(
             {      v0 * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5),
                    v_bar * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5),
                    rho * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5),
                    kappa * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5),
                    sigma * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5)
            }, (((double) rand() / (RAND_MAX))),0.02);
        double x=(((double) rand() / (RAND_MAX)) - 0.5)*500;
        std::vector<std::complex<double>> correct =GetCuiGradient(test.p,-x,test.tau);
        std::vector<std::complex<double>> wrong = test.chf_grad(x);
        double error=0;
        for(int j=0;j<5;j++) error+=std::abs(correct[0]-wrong[0]);
        std::cout<<"v_0: "<<test.p.v_0<<"\tv_m: "<<test.p.v_m<<"\trho: "<<test.p.rho<<"\tkappa: "<<test.p.kappa<<"\tsigma: "<<test.p.sigma<<"\ttau: "<<test.tau<<"\tx: "<<x<<"\tdiff: "<<std::sqrt(error)<<'\n';
        assert(std::sqrt(error)<DISTR_TOLERANCE);
    }
    std::cout<<"Distribution Test passed"<<std::endl;
}
//#################################################################Model-Calibration#########################################################################################################
std::vector<double> expiries = {
    0.119047619047619,
    0.238095238095238,
    0.357142857142857,
    0.476190476190476,
    0.595238095238095,
    0.714285714285714,
    1.07142857142857,
    1.42857142857143};
std::vector<std::vector<double>> strikes = {
    {0.9371, 0.9956, 1.0427, 1.2287, 1.3939},
    {0.8603, 0.9868, 1.0463, 1.2399, 1.4102},
    {0.8112, 0.9728, 1.0499, 1.2485, 1.4291},
    {0.7760, 0.9588, 1.0530, 1.2659, 1.4456},
    {0.7470, 0.9464, 1.0562, 1.2646, 1.4603},
    {0.7216, 0.9358, 1.0593, 1.2715, 1.4736},
    {0.6699, 0.9175, 1.0663, 1.2859, 1.5005},
    {0.6137, 0.9025, 1.0766, 1.3046, 1.5328}
};
std::vector<swift_parameters> params={{5, 32, 5.6568542494923806, -1.8343994686679572, 1.5720210785285174, -58, 50, 128},
                                        {5, 32, 5.6568542494923806, -2.6596460616572828, 2.4759124462331239, -85, 79, 256},
                                        {5, 32, 5.6568542494923806, -3.3440060583236764, 3.2104875432489965, -107, 102, 256},
                                        {5, 32, 5.6568542494923806, -3.9452944230550919, 3.8494203406057474, -126, 123, 256},
                                        {5, 32, 5.6568542494923806, -4.4898573506150337, 4.4267150742095884, -143, 141, 256},
                                        {5, 32, 5.6568542494923806, -4.9920925407415853, 4.9592398930278314, -159, 158, 256},
                                        {5, 32, 5.6568542494923806, -6.3157233304280682, 6.353408918854841, -202, 203, 512},
                                        {5, 32, 5.6568542494923806, -7.4606624599663149, 7.5789582584617987, -238, 242, 512}};
double p[5]={
    .09,          // v0
    .09,          // v_bar
    .04,          // rho
    1.,            // kappa
    1.            // sigma
};
std::unique_ptr<adata_s> get_adata(const double S, const std::vector<double>& expiries_arg,const std::vector<std::vector<double>>& strikes_arg, const HParams heston_params,bool predefined){
    std::unique_ptr<adata_s> adata=std::unique_ptr<adata_s>(new adata_s{S,nullptr,*(new std::list<expiry_data>())});
    for (unsigned int i = 0; i < expiries_arg.size(); ++i){
        std::vector<ffloat> cur=strikes_arg[i];
        options_chain * opt_chain=new options_chain(static_cast<unsigned int>(expiries_arg[i])/trading_days,expiries_arg[i]);
        for(auto c: cur){
            option * new_opt =new option();
            new_opt->volume=1;
            new_opt->strike=c;
            new_opt->price=1.;
            opt_chain->options->push_back(*new_opt);
        }
        opt_chain->min_strike=cur[0];
        opt_chain->max_strike=cur[cur.size()-1];
        HDistribution *new_distr=new HDistribution(heston_params,expiries_arg[i],0.02);
        swift_parameters dynamic_parameters(*new_distr,adata->S,*opt_chain);
        swift_parameters& new_swift_parameters=predefined ? params[i] : dynamic_parameters;
        SWIFT* pricing_method=new SWIFT(new_swift_parameters);//std::shared_ptr(current);
        adata->exp_list.emplace_back(*opt_chain,new_distr,pricing_method);
    }
    return adata;
}
#define PRICING_TOLERANCE 1e-9
void pricing_test(){
    double diff=0.;
    std::vector<ffloat> prices={
        0.079676812094469612,
        0.042586263756033402,
        0.023642876097251266,
        0.0019447635553313004,
        0.00022457675334788794,
        0.15394308829999556,
        0.062423699148410512,
        0.03716156027004304,
        0.0075329906749080954,
        0.0021671474663124877,
        0.20438572427747337,
        0.081424876885654141,
        0.047007477654992851,
        0.013394276263081869,
        0.0051665880568379507,
        0.24157958905122354,
        0.09899570130848076,
        0.055258126136188711,
        0.017997398670040403,
        0.0084472770368514277,
        0.27251527545827303,
        0.11480331435058673,
        0.062447187930262341,
        0.023683757123971076,
        0.011714319060476305,
        0.29968779707061638,
        0.12892214079993489,
        0.068994812796204188,
        0.02829225908895841,
        0.01488391596622375,
        0.3569911166621863,
        0.16022472995208151,
        0.086978384558870039,
        0.041295845512103642,
        0.024161759165286907,
        0.41608092405221836,
        0.18693437912935496,
        0.10163212168506526,
        0.052461532814837883,
        0.032268432242168528};
    std::unique_ptr<adata_s> adata_ptr=get_adata(1.,expiries, strikes,{p[0],p[1],p[2],p[3],p[4]}, true);
    ffloat * x=(ffloat*) malloc(sizeof(ffloat)*40);
    get_prices_for_levmar(&p[0], x, 5, 40, (void *) &(*adata_ptr));
    for(int i=0;i<40;i++){
        double error=std::fabs(x[i]-prices[i]);
        std::cout<<"x: "<<x[i]<<"\tp: "<<prices[i]<<"\tdiff: "<<error<<'\n';
        diff+=error;
    }
    assert(diff<PRICING_TOLERANCE);
    std::cout<<"Pricing test passed"<<std::endl;
}
#define GRADIENT_TOLERANCE 1e-9
void gradient_test(){
    double diff=0.;
        std::vector<ffloat> grad={6.0887682383150116e-05,
                                0.01009182616033304,
                                -0.001422153712575771,
                                -0.0052144883688743328,
                                0.1622664757551231,
                                0.00019059712892998947,
                                0.014367167944521392,
                                -0.0038230359319571801,
                                -0.00026348078537531157,
                                0.22358069544554257,
                                0.0001243331501459066,
                                0.012962315883333041,
                                -0.0026336197669243888,
                                0.0047260314836243779,
                                0.20557573274688637,
                                -5.9753433147245287e-05,
                                0.0021031499469935807,
                                0.0013282359979569184,
                                0.0032813563925235176,
                                0.038629111670166423,
                                -1.7507894271232844e-05,
                                0.0003026256761754147,
                                0.00042754784119065225,
                                0.0007090133895697175,
                                0.0060337930478490363,
                                -0.00010109555807403782,
                                0.015522478525322046,
                                0.0005310649843597307,
                                -0.0084268534027068034,
                                0.12957038021122541,
                                0.00084409938092525505,
                                0.040269636611385821,
                                -0.0081020288753043384,
                                -0.0014793027431618986,
                                0.29552667707773,
                                0.00068978146591164595,
                                0.038407499235951138,
                                -0.0068864772208932614,
                                0.0072585866071890204,
                                0.28686172501536633,
                                -0.00026252489005407963,
                                0.01171971690728809,
                                0.0024816357688686153,
                                0.0094110571536697649,
                                0.10463610152052862,
                                -0.00020190950059931933,
                                0.0038173484133046001,
                                0.0023133200997451833,
                                0.0045947094374290728,
                                0.03739938463348761,
                                -0.00029652546401583504,
                                0.021186983600339795,
                                0.0013750309095158131,
                                -0.009717046739688982,
                                0.11964948927773626,
                                0.0017176406300190663,
                                0.069778508989456381,
                                -0.011017963220314564,
                                -0.0036722895906964766,
                                0.32850208612571719,
                                0.0016297018402991529,
                                0.070431184675385089,
                                -0.010639071004535499,
                                0.0091294850313919634,
                                0.33483916920474344,
                                -0.00040662733381871525,
                                0.027413672689554153,
                                0.0020127421159107377,
                                0.01460550625271187,
                                0.15648310653349168,
                                -0.00053345940166278616,
                                0.011289443288832027,
                                0.003774241733427627,
                                0.0092447362384140905,
                                0.07204037129236246,
                                -0.00048487568237889195,
                                0.027703202063868119,
                                0.0017494596153002527,
                                -0.010648728991968443,
                                0.11706553370902192,
                                0.0026128883245209004,
                                0.099454606059466391,
                                -0.012856656432388095,
                                -0.0057557658913496602,
                                0.34179817772695137,
                                0.0028139209661547661,
                                0.10622666517857729,
                                -0.013811917234721213,
                                0.010634906596537973,
                                0.36515595503801829,
                                -0.00044001435338978132,
                                0.045945304640325778,
                                0.0011085103980652492,
                                0.018714649018874398,
                                0.1893609171102367,
                                -0.00089853034744459786,
                                0.022202823862855086,
                                0.0044337545861627961,
                                0.013699960722508961,
                                0.10277861910442491,
                                -0.00066258982660075583,
                                0.034406191436722709,
                                0.0019591262143772506,
                                -0.011332692849232053,
                                0.11515126420375622,
                                0.0034835048950391592,
                                0.12857759040110306,
                                -0.014118994531699089,
                                -0.0074094651083295273,
                                0.34591456716843927,
                                0.0041387978240080107,
                                0.14403121305187261,
                                -0.01649439149813477,
                                0.011996725869116301,
                                0.38442080038735577,
                                -0.00015309583873739715,
                                0.0718326535827878,
                                -0.00081333310763814434,
                                0.02258212525054467,
                                0.22519701702181644,
                                -0.0012105403856448283,
                                0.036179416807485924,
                                0.0044878946805707743,
                                0.017798629739292332,
                                0.12880052004981041,
                                -0.00082961389580639389,
                                0.040983691938186143,
                                0.0021021053915233324,
                                -0.01182763131981346,
                                0.1127729999934402,
                                0.0043326060861786893,
                                0.15718926372231473,
                                -0.015096475918122375,
                                -0.008642130678579784,
                                0.34572653078530097,
                                0.0055438716932617125,
                                0.18288390299297339,
                                -0.018816401421951172,
                                0.013254095333921422,
                                0.39656435329702744,
                                0.00025765617248749379,
                                0.098897874554757584,
                                -0.0024838269394367031,
                                0.025838401458442746,
                                0.2483016099909742,
                                -0.0014296651834192595,
                                0.052809101906874913,
                                0.0041697971039235922,
                                0.021586277717352213,
                                0.15045329378680108,
                                -0.0011704248631786412,
                                0.063828470352787756,
                                0.0020836161253029936,
                                -0.013107116819623312,
                                0.11026888983742654,
                                0.0070764368715066249,
                                0.24522986441583147,
                                -0.017965550011084831,
                                -0.010092208052498946,
                                0.34101797638216724,
                                0.0099059732604421247,
                                0.3005757492296996,
                                -0.024384078320241056,
                                0.016321992971710341,
                                0.40929106415535943,
                                0.0024113239826787924,
                                0.19314682141287215,
                                -0.0076072155732986894,
                                0.033961547807685462,
                                0.29291928335576756,
                                -0.0013788932223915112,
                                0.11676103917807525,
                                0.0020276924471335309,
                                0.031715203627352598,
                                0.19784895628220139,
                                -0.0015487506643329322,
                                0.077086883579066387,
                                0.0024510287956627554,
                                -0.013160026523107609,
                                0.095601550969006507,
                                0.0095488642081517586,
                                0.32678279093067242,
                                -0.019924378831838199,
                                -0.010806406891092421,
                                0.32662447748208967,
                                0.013995824693364179,
                                0.41323259877593316,
                                -0.028396705182887921,
                                0.019694514760990905,
                                0.4038676861332352,
                                0.0051168734628484151,
                                0.29266055294591559,
                                -0.011818563276080525,
                                0.040641123734601339,
                                0.31052466854997196,
                                -0.00059682554659867284,
                                0.19035093225791697,
                                -0.00030980371835320395,
                                0.040232991901821259,
                                0.22181984888582221};
    std::vector<ffloat> index_map={4,1,3,0,2};
    std::unique_ptr<adata_s> adata_ptr=get_adata(1.,expiries, strikes,{p[0],p[1],p[2],p[3],p[4]},true);
    ffloat * jac=(ffloat*) malloc(sizeof(ffloat)*200);
    get_jacobian_for_levmar(&p[0], jac, 5, 40, (void *) &(*adata_ptr));
    for(int i=0;i<40;i++){
        double error=0.;
        for(int j=0;j<5;j++) error+=std::fabs(grad[5*i+index_map[j]]-jac[5*i+j]);
        std::cout<<"dv0(ER): "<<grad[5*i+4]<<"\tdv0(ME): "<<jac[5*i]
        <<"\tdvm(ER): "<<grad[5*i+1]<<"\tdvm(ME): "<<jac[5*i+1]
        <<"\tdrho(ER): "<<grad[5*i+3]<<"\tdrho(ME): "<<jac[5*i+2]
        <<"\tdkappa(ER): "<<grad[5*i]<<"\tdkappa(ME): "<<jac[5*i+3]
        <<"\tdsigma(ER): "<<grad[5*i+2]<<"\tdsigma(ME): "<<jac[5*i+4]
        <<"\ndiff: "<<error<<'\n';
        diff+=error;
    }
    assert(diff<GRADIENT_TOLERANCE);
    std::cout<<"Gradient test passed"<<std::endl;
}
#define LEVMAR_TOLERANCE 1e-2
void levmar_test(){
    double p2[5];
    std::unique_ptr<adata_s> adata_ptr=get_adata(1.,expiries, strikes,{p[0],p[1],p[2],p[3],p[4]},false);
    ffloat * x=(ffloat*) malloc(sizeof(ffloat)*40);
    get_prices_for_levmar(&p[0], x, 5, 40, (void *) &(*adata_ptr));
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU;
    // stopping thresholds for
    opts[1]=1e-10;       // ||J^T e||_inf
    opts[2]=1e-10;       // ||Dp||_2
    opts[3]=1e-10;       // ||e||_2
    opts[4]= LM_DIFF_DELTA; // finite difference if used
    int retval;
    p2[0]=p[0]+.2;p2[1]=p[1]+.2,p2[2]=p[2]+.2;p2[3]=p[3]+.2;p2[4]=p[4]+.2;
    adata_ptr=get_adata(1.,expiries, strikes,{p2[0],p2[1],p2[2],p2[3],p2[4]},false);
    retval=dlevmar_der(get_prices_for_levmar, get_jacobian_for_levmar, p2, x, 5, 40, 100, opts, info, NULL, NULL, (void*) &(*adata_ptr));
    double diff=0.;
    for(int i=0;i<5; i++) diff+=std::fabs(p[i]-p2[i]);
    std::cout<<"# iter: "<<retval<<"\tv_0: "<<p2[0]<<"\tv_m: "<<p2[1]<<"\trho: "<<p2[2]<<"\tkappa: "<<p2[3]<<"\tsigma: "<<p2[4]<<"\tinital error: "<<info[0]<<"\te: "<<info[1]<<"\tparams diff: "<<diff<<"\treason: "<<info[6]<<'\n';
    assert(diff<LEVMAR_TOLERANCE &&"Levmar test passed");
    std::cout<<"Levmar test passed"<<std::endl;
}
#define SAMPLE_SIZE (1<<28)
#define BLOCK_SIZE (1<<20)
#define RNG_TOLERANCE 1e-4
#define PROG_BAR_WIDTH 100
#define PROG_CHAR '#'
void update_prog(char* prog_bar, unsigned int i){
    double prog= ((double)i)/((double)SAMPLE_SIZE);
    int visible_prog=(int)(prog*PROG_BAR_WIDTH);
    printf("\r[%.*s%*s] %3d%%", visible_prog, prog_bar, PROG_BAR_WIDTH-visible_prog, " ",(int) (prog*100));
    fflush(stdout);
}
void rng_test(){
    RNG shishua(BLOCK_SIZE, 1);
    ffloat mean_uniform=0.;
    ffloat emp_uniform=0.;
    ffloat mean_gaussian=0.;
    ffloat emp_gaussian=0.;
    std::cout<<"Start RNG test, sample size: "<<SAMPLE_SIZE<<", block size: "<<BLOCK_SIZE<<std::endl;
    //ffloat* uniform_array=(ffloat*) malloc(sizeof(ffloat)*SAMPLE_SIZE);
    //ffloat* gaussian_array=(ffloat*) malloc(sizeof(ffloat)*SAMPLE_SIZE);
    char prog_bar[PROG_BAR_WIDTH+1];
    for(unsigned int i=0;i<PROG_BAR_WIDTH;i++) prog_bar[i]=PROG_CHAR;
    prog_bar[PROG_BAR_WIDTH]='\0';
    //std::cout<<"prog bar initialized"<<prog_bar<<std::endl;
    for(unsigned int i=0;i<SAMPLE_SIZE;i++){
        mean_uniform+=shishua.get_urand();
        mean_gaussian+=shishua.get_grand();
        if(!(i&((SAMPLE_SIZE>>7)-1))) update_prog(prog_bar,i);
    }
    printf("\r[%s] 100%%\n",prog_bar);
    mean_uniform/=SAMPLE_SIZE;
    mean_gaussian/=SAMPLE_SIZE;
    RNG shishua2(BLOCK_SIZE, 1); //reinitialize shishua since RNG is deterministic we can do this and will get the exact same (pseudo) random number
    for(unsigned int i=0;i<SAMPLE_SIZE;i++){
        double uni_rand=shishua2.get_urand();
        double gaussian_rand=shishua2.get_grand();
        emp_uniform+=(uni_rand-mean_uniform)*(uni_rand-mean_uniform);
        emp_gaussian+=(gaussian_rand-mean_gaussian)*(gaussian_rand-mean_gaussian);
        if(!(i&((SAMPLE_SIZE>>7)-1))) update_prog(prog_bar,i);
    }
    printf("\r[%s] 100%%\n",prog_bar);
    emp_uniform/=(SAMPLE_SIZE-1);
    emp_gaussian/=(SAMPLE_SIZE-1);
    std::cout<<"emp_uniform: "<<emp_uniform<<"\tmean_uniform: "<<mean_uniform<<"\temp_gaussian: "<<emp_gaussian<<"\tmean_gaussian: "<<mean_gaussian<<std::endl;
    std::cout<<"emp_uniform_error: "<<std::fabs(emp_uniform-1./12.)<<"\tmean_uniform_error: "<<std::fabs(mean_uniform-.5)<<"\temp_gaussian_error: "<<std::fabs(emp_gaussian-1.)<<"\tmean_gaussian_error: "<<std::fabs(mean_gaussian-0.)<<std::endl;
    assert(std::fabs(emp_gaussian-1.)<RNG_TOLERANCE&&std::fabs(mean_gaussian-0.)<RNG_TOLERANCE&&
            std::fabs(mean_uniform-.5)<RNG_TOLERANCE&&std::fabs(emp_uniform-1./12.)<RNG_TOLERANCE);
}
void simulation_test(){
    std::vector<ffloat> deltas={1/32,1/16,1/8,1/4,1/2,1};
    std::vector<ffloat> strikes={70.,100.,140.};
    std::vector<HParams> ps={{0.04,0.04,-0.9, 0.4,1.},
                             {0.04,0.04,-0.5, 0.3,1.},
                             {0.09,0.09,-0.3, 1. ,1.}};
    std::vector<ffloat> years={5.,10.,15.};
    for(int i=0;i<3;i++){
        std::cout<<"batch i: "<<i<<std::endl;
        std::list<options_chain>all_chains;
        all_chains.emplace_back(years[i]/trading_days,years[i]);
        options_chain opts=*all_chains.begin();
        for(ffloat strike: strikes) opts.options->push_back({0.,0.,strike,0});
        for(ffloat delta: deltas){
            HSimulation::PricingTool<ffloat,EuropeanCallNonAdaptive>my_pricing_tool(1);
            my_pricing_tool.price(ps[i],100,all_chains,1e+6,3,years[i]/delta);
        }
    }
}

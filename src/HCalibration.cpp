/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"HCalibration.h"
#include"BSM.h"
#include"SWIFT.h"
#include<iostream>
#include<levmar/levmar.h>
/*typedef struct {
    ffloat mu;
    ffloat theta;
    ffloat kappa;
    ffloat xi;
    ffloat rho;
} HestonParams;
*/
typedef std::numeric_limits< double > dbl;

typedef struct ED:traced<ED>{
    const options_chain& opts;
    HDistribution* distr;
    SWIFT* pricing_method;
    //unsigned int underpriced;
    ED(const options_chain& opts,HDistribution* init_distr, SWIFT* init_pricing_method): opts(opts),distr(init_distr),pricing_method(init_pricing_method){}
    //ED(const options_chain& opts,HParams p, ffloat expi): opts(opts),distr(new HDistribution(p,expi)){}
    ~ED(){delete distr; delete pricing_method;}
} expiry_data;
typedef struct AS{
    ffloat S;
    ffloat* real_prices;
    std::list<expiry_data>& exp_list;
    ~AS(){delete &exp_list;} //moved when push_back 
} adata_s;
#define NEWOLD_METHOD_RATIO 0.1
void update_adata(ffloat *p, adata_s * adata);

void get_prices_for_levmar(ffloat *p, ffloat *x, int m, int n_observations, void * adata);
void get_jacobian_for_levmar(ffloat *p, ffloat *jac, int m, int n_observations, void * adata);

void update_adata(ffloat *p, adata_s * adata){
    const HParams new_params={p[0],p[1],p[2],p[3],p[4]};
    //std::shared_ptr<SWIFT> current=nullptr; delete
    for(auto &exp_data : adata->exp_list){
        if(!(exp_data.distr->p==new_params)){
            exp_data.pricing_method->flush_cache();
            //std::cout<<"delete distr\n";
            exp_data.distr->p=new_params;//}
            //std::cout<<"new params at: "<<&(exp_data.distr->p)<<"\tv_0: "<<exp_data.distr->p.v_0<<"\tv_m: "<<exp_data.distr->p.v_m<<"\trho: "<<exp_data.distr->p.rho<<"\tkappa"<<exp_data.distr->p.kappa<<"\tsigma: "<<exp_data.distr->p.sigma<<'\n';
        //std::cout<<"cur params at: "<<&(exp_data.distr->p)<<"\tv_0: "<<(float)exp_data.distr->p.v_0<<"\tv_m: "<<exp_data.distr->p.v_m<<"\trho: "<<exp_data.distr->p.rho<<"\tkappa"<<exp_data.distr->p.kappa<<"\tsigma: "<<exp_data.distr->p.sigma<<'\n';
        //std::cout<<"new SWIFT\n";
        auto new_swift_parameters=SWIFT::get_parameters(*exp_data.distr,adata->S,exp_data.opts);
        //if (new_swift_parameters->m>exp_data.pricing_method->my_params.m 
            //||new_swift_parameters->J<NEWOLD_METHOD_RATIO*exp_data.pricing_method->my_params.J
        //){
            //delete exp_data.pricing_method;
        //    if(current==nullptr||new_swift_parameters->m>current->my_params.m ||new_swift_parameters->J<NEWOLD_METHOD_RATIO*current->my_params.J)
        //        current=std::make_shared<SWIFT>(*new_swift_parameters);
            //std::cout<<"old m: "<<new_swift_parameters->m<<"\tnew m: "<<exp_data.pricing_method->my_params.m<<'\n'; 
            delete exp_data.pricing_method;
            exp_data.pricing_method=new SWIFT(*new_swift_parameters);//std::shared_ptr(current);
        }
    }
}
void get_jacobian_for_levmar(ffloat *p, ffloat *jac, int m, int n_observations, void * adata){
    std::cout<<"get gradient\tv_0: "<<p[0]<<"\tv_m: "<<p[1]<<"\trho: "<<p[2]<<"\tkappa: "<<p[3]<<"\tsigma: "<<p[4]<<'\n';
    adata_s * my_adata=static_cast<adata_s*>(adata);
    update_adata(p,my_adata);
    ffloat* jac2=jac;
    for(auto &exp_data : my_adata->exp_list) exp_data.pricing_method->price_opts_grad(*exp_data.distr,my_adata->S,exp_data.opts, &jac2,jac+n_observations*m);
    if(jac2<jac+n_observations*m) throw std::runtime_error("Gradient buffer too large");
}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void get_prices_for_levmar(ffloat *p, ffloat *x, int m, int n_observations, void * adata){
    std::cout<<"get prices\tv_0: "<<p[0]<<"\tv_m: "<<p[1]<<"\trho: "<<p[2]<<"\tkappa: "<<p[3]<<"\tsigma: "<<p[4]<<'\n';
    adata_s * my_adata=static_cast<adata_s*>(adata);
    update_adata(p,my_adata);
    ffloat* x2=x;
    for(auto &exp_data : my_adata->exp_list) exp_data.pricing_method->price_opts(*exp_data.distr,my_adata->S,exp_data.opts, &x2,x+n_observations);
    if(x2<x+n_observations) throw std::runtime_error("Pricing buffer too large");
    //unsigned int underpriced=0;
    //if(my_adata->real_prices) for(x2=x;x2<x+n_observations;x2++) if(*x2<my_adata->real_prices[x2-x]) underpriced++;
    //std::cout<<"share of underpriced: "<<static_cast<ffloat>(underpriced)/n_observations<<'\n';
}
#pragma GCC diagnostic pop

std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<options_chain>& market_data){
    adata_s adata={S,nullptr,*(new std::list<expiry_data>())};
    //HDistribution current_distribution;
    //std::shared_ptr<SWIFT> current;
    ffloat avg_iv=avg_imp_vol(S, market_data);
    // We start in a Black-Scholes Model see p.106 lecture notes
    ffloat yearly_avg_iv = avg_iv*avg_iv;
    //std::cout<<"yearly_avg_iv: "<<yearly_avg_iv<<'\n';
    ffloat p[5];
    p[0]=yearly_avg_iv;
    p[1]=yearly_avg_iv;
    p[2]=0.04;
    p[3]=1.;
    p[4]=1.;
    unsigned int n_observations_cur=0;
    std::cout<<"start levmar setup\n";
    for(const auto &opts: market_data){
        //if(opts->time_to_expiry<=EXP_LB) continue;
        if(opts.options->size()>0){
            n_observations_cur+=opts.options->size();
            HDistribution *new_distr=new HDistribution({p[0],p[1],p[2],p[3],p[4]},opts.time_to_expiry,0.0005);
            auto new_swift_parameters=SWIFT::get_parameters(*new_distr,S,opts);
            //if (current==nullptr|| new_swift_parameters->m>current->my_params.m || new_swift_parameters->J<NEWOLD_METHOD_RATIO*current->my_params.J){
                //std::cout<<"is here the free1?\n";
                //current=std::make_shared<SWIFT>(*new_swift_parameters);
                //std::cout<<"is here the free?\n";
            //}
            SWIFT* pricing_method=new SWIFT(*new_swift_parameters);//std::shared_ptr(current);
            adata.exp_list.emplace_back(opts,new_distr,pricing_method);
        }
    }
    //std::cout<<"allocate x\n";
    ffloat * x=(ffloat*) malloc(sizeof(ffloat)*(n_observations_cur+1));
    ffloat * x2=x;
    for(std::list<options_chain>::const_iterator opts = market_data.begin(); opts != market_data.end(); opts++) for(auto e :*opts->options) *(x2++)=e.price;
    adata.real_prices=x;
    std::cout<<"setup completed\t# observations: "<<n_observations_cur<<'\n';
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU;
    // stopping thresholds for
    opts[1]=1E-10;       // ||J^T e||_inf
    opts[2]=1E-10;       // ||Dp||_2
    opts[3]=1E-10;       // ||e||_2
    opts[4]= LM_DIFF_DELTA; // finite difference if used
    std::vector<std::string> msg={{"1 -stopped by small gradient J^T e"},
                                {"2 - stopped by small Dp"},
                                {"3 - stopped by itmax"},
                                {"4 - singular matrix. Restart from current p with increased \\mu"},
                                {"5 - no further error reduction is possible. Restart with increased mu"},
                                {"6 - stopped by small ||e||_2"},
                                {"7 - stopped by invalid (i.e. NaN or Inf) \"func\" values; a user error"}};
    int iter=dlevmar_der(get_prices_for_levmar, get_jacobian_for_levmar, p, x, 5, n_observations_cur, 100, opts, info, NULL, NULL, (void*) &adata);
    if(iter<0) throw std::runtime_error("levmar failed!");
    auto to_calib=std::unique_ptr<HParams>(new HParams({p[0],p[1],p[2],p[3],p[4]}));
    std::cout<<"# iter: "<<iter<<"\tv_0: "<<to_calib->v_0<<"\tv_m: "<<to_calib->v_m<<"\trho: "<<to_calib->rho<<"\tkappa: "<<to_calib->kappa<<"\tsigma: "<<to_calib->sigma<<"\tinital e: "<<info[0]<<"\te: "<<info[1]<<"\treason: "<<msg[info[6]-1]<<'\n';
    if(info[6]!=6.&&info[6]!=2) throw std::runtime_error("levmar failed! ");
    return to_calib;
}

void pricing_test(){
    //yearly_risk_free=0.02;
    std::cout.precision(dbl::max_digits10);
    adata_s adata={1.,nullptr,*(new std::list<expiry_data>())};
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
    std::vector<double> expiries = {
    0.119047619047619,
    0.238095238095238,
    0.357142857142857,
    0.476190476190476,
    0.595238095238095,
    0.714285714285714,
    1.07142857142857,
    1.42857142857143};
    std::vector<std::vector<double>> K_over_S = {
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
    double kappa = 1;           // |  mean reversion rate
    double v_bar = 0.09;          // |  long term variance
    double sigma = 1;          // |  variance of volatility
    double rho = 0.04;            // |  correlation between spot and volatility
    double v0 = 0.09;
    double p[5];
    p[0]=v0;p[1]=v_bar,p[2]=rho;p[3]=kappa;p[4]=sigma;
    for (unsigned int index = 0; index < expiries.size(); ++index)
    {
        std::vector<ffloat> cur=K_over_S[index];
        options_chain * opt_chain=new options_chain(static_cast<unsigned int>(expiries[index])/trading_days,expiries[index]);
        for(auto c: cur){
            option * new_opt =new option();
            new_opt->volume=1;
            new_opt->strike=c;
            new_opt->price=1.;
            opt_chain->options->push_back(*new_opt);
        }
        opt_chain->min_strike=cur[0];
        opt_chain->max_strike=cur[cur.size()-1];
        HDistribution *new_distr=new HDistribution({v0,v_bar,rho,kappa,sigma},expiries[index],0.02);
        //auto new_swift_parameters=SWIFT::get_parameters(*new_distr,adata.S,*opt_chain);
        SWIFT* pricing_method=new SWIFT(params[index]);//std::shared_ptr(current);
        adata.exp_list.emplace_back(*opt_chain,new_distr,pricing_method);
    }
    ffloat * x=(ffloat*) malloc(sizeof(ffloat)*40);
    get_prices_for_levmar(&p[0], x, 5, 40, (void *) &adata);
    for(int i=0;i<40;i++){
        std::cout<<"x: "<<x[i]<<"\tp: "<<prices[i]<<"\tdiff: "<<std::fabs(x[i]-prices[i])<<'\n';
    }
}
void gradient_test(){
    //yearly_risk_free=0.02;
    adata_s adata={1.,nullptr,*(new std::list<expiry_data>())};
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
    std::vector<double> expiries = {
    0.119047619047619,
    0.238095238095238,
    0.357142857142857,
    0.476190476190476,
    0.595238095238095,
    0.714285714285714,
    1.07142857142857,
    1.42857142857143};
    std::vector<std::vector<double>> K_over_S = {
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
    double kappa = 1;           // |  mean reversion rate
    double v_bar = 0.09;          // |  long term variance
    double sigma = 1;          // |  variance of volatility
    double rho = 0.04;            // |  correlation between spot and volatility
    double v0 = 0.09;
    double p[5];
    p[0]=v0;p[1]=v_bar,p[2]=rho;p[3]=kappa;p[4]=sigma;
    for (unsigned int index = 0; index < expiries.size(); ++index)
    {
        std::vector<ffloat> cur=K_over_S[index];
        options_chain * opt_chain=new options_chain(static_cast<unsigned int>(expiries[index])/trading_days,expiries[index]);
        for(auto c: cur){
            option * new_opt =new option();
            new_opt->volume=1;
            new_opt->strike=c;
            new_opt->price=1.;
            opt_chain->options->push_back(*new_opt);
        }
        opt_chain->min_strike=cur[0];
        opt_chain->max_strike=cur[cur.size()-1];
        HDistribution *new_distr=new HDistribution({v0,v_bar,rho,kappa,sigma},expiries[index],0.02);
        //auto new_swift_parameters=SWIFT::get_parameters(*new_distr,adata.S,*opt_chain);
        SWIFT* pricing_method=new SWIFT(params[index]);//std::shared_ptr(current);
        adata.exp_list.emplace_back(*opt_chain,new_distr,pricing_method);
    }
    ffloat * jac=(ffloat*) malloc(sizeof(ffloat)*200);
    get_jacobian_for_levmar(&p[0], jac, 5, 40, (void *) &adata);
    for(int i=0;i<40;i++){
        std::cout<<"dv0(ER): "<<grad[5*i+4]<<"\tdv0(ME): "<<jac[5*i]
        <<"\tdvm(ER): "<<grad[5*i+1]<<"\tdvm(ME): "<<jac[5*i+1]
        <<"\tdrho(ER): "<<grad[5*i+3]<<"\tdrho(ME): "<<jac[5*i+2]
        <<"\tdkappa(ER): "<<grad[5*i]<<"\tdkappa(ME): "<<jac[5*i+3]
        <<"\tdsigma(ER): "<<grad[5*i+2]<<"\tdsigma(ME): "<<jac[5*i+4]
        <<"\ndiff: "<<std::fabs(grad[5*i+4]-jac[5*i])+
                    std::fabs(grad[5*i+1]-jac[5*i+1])+
                    std::fabs(grad[5*i+3]-jac[5*i+2])+
                    std::fabs(grad[5*i]-jac[5*i+3])+
                    std::fabs(grad[5*i+2]-jac[5*i+4])<<'\n';
    }
}
void levmar_test(){
    //yearly_risk_free=0.02;
    adata_s adata={1.,nullptr,*(new std::list<expiry_data>())};
    std::vector<double> expiries = {
    0.119047619047619,
    0.238095238095238,
    0.357142857142857,
    0.476190476190476,
    0.595238095238095,
    0.714285714285714,
    1.07142857142857,
    1.42857142857143};
    std::vector<std::vector<double>> K_over_S = {
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
    double kappa = 1;           // |  mean reversion rate
    double v_bar = 0.09;          // |  long term variance
    double sigma = 1;          // |  variance of volatility
    double rho = 0.04;            // |  correlation between spot and volatility
    double v0 = 0.09;
    double p[5];double p2[5];
    p[0]=v0;p[1]=v_bar,p[2]=rho;p[3]=kappa;p[4]=sigma;
    p2[0]=v0+.2;p2[1]=v_bar+.2,p2[2]=rho+.2;p2[3]=kappa+.2;p2[4]=sigma+.2;
    for (unsigned int index = 0; index < expiries.size(); ++index)
    {
        std::vector<ffloat> cur=K_over_S[index];
        options_chain * opt_chain=new options_chain(static_cast<unsigned int>(expiries[index])/trading_days,expiries[index]);
        for(auto c: cur){
            option * new_opt =new option();
            new_opt->volume=1;
            new_opt->strike=c;
            new_opt->price=1.;
            opt_chain->options->push_back(*new_opt);
        }
        opt_chain->min_strike=cur[0];
        opt_chain->max_strike=cur[cur.size()-1];
        HDistribution *new_distr=new HDistribution({v0,v_bar,rho,kappa,sigma},expiries[index],0.02);
        auto new_swift_parameters=SWIFT::get_parameters(*new_distr,adata.S,*opt_chain);
        SWIFT* pricing_method=new SWIFT(*new_swift_parameters);//std::shared_ptr(current);
        adata.exp_list.emplace_back(*opt_chain,new_distr,pricing_method);
    }
    ffloat * x=(ffloat*) malloc(sizeof(ffloat)*40);
    get_prices_for_levmar(&p[0], x, 5, 40, (void *) &adata);
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU;
    // stopping thresholds for
    opts[1]=1E-10;       // ||J^T e||_inf
    opts[2]=1E-10;       // ||Dp||_2
    opts[3]=1E-10;       // ||e||_2
    opts[4]= LM_DIFF_DELTA; // finite difference if used
    int retval;
    adata.exp_list.clear();
    for (unsigned int index = 0; index < expiries.size(); ++index)
    {
        std::vector<ffloat> cur=K_over_S[index];
        options_chain * opt_chain=new options_chain(static_cast<unsigned int>(expiries[index])/trading_days,expiries[index]);
        for(auto c: cur){
            option * new_opt =new option();
            new_opt->volume=1;
            new_opt->strike=c;
            new_opt->price=1.;
            opt_chain->options->push_back(*new_opt);
        }
        opt_chain->min_strike=cur[0];
        opt_chain->max_strike=cur[cur.size()-1];
        HDistribution *new_distr=new HDistribution({p2[0],p2[1],p2[2],p2[3],p2[4]},expiries[index],0.02);
        auto new_swift_parameters=SWIFT::get_parameters(*new_distr,adata.S,*opt_chain);
        SWIFT* pricing_method=new SWIFT(*new_swift_parameters);//std::shared_ptr(current);
        adata.exp_list.emplace_back(*opt_chain,new_distr,pricing_method);
    }
    retval=dlevmar_der(get_prices_for_levmar, get_jacobian_for_levmar, p2, x, 5, 40, 100, opts, info, NULL, NULL, (void*) &adata);
    std::cout<<"# iter: "<<retval<<"\tv_0: "<<p2[0]<<"\tv_m: "<<p2[1]<<"\trho: "<<p2[2]<<"\tkappa: "<<p2[3]<<"\tsigma: "<<p2[4]<<"\tinital error: "<<info[0]<<"\te: "<<info[1]<<"\treason: "<<info[6]<<'\n';
}

#include "asano_model.h"

#include <SED_setting_util.h>

// AGBdust:Zhukovska(2008), SNdust:Bianchi&Schneider(2007) or Nozawa et al.(2007)
// AGBmetal:Hoeck&groenewegen(1997), SNmetal:Woosley&Weaver(1995)
// AGBremnant:Hoeck&groenewegen(1997), SNremnant:Woosley&Weaver(1995)
using SEDSetting = my_util::SED_util::SEDSettingUtil;

int main() {
    // 使われてないっぽい
    // const auto is_infall      = true ;
    // const auto M_total_gal    = 1e11;
    // const auto tau_SF         = 3e9;
    // const auto schmidt_index  = 1.0;
    // const auto M_total_infall = M_total_gal;
    // const auto tau_infall     = 6e9;    //
    // const auto age_max        = 13e9;
    //

    const auto asano_model = asano_model::AsanoModel();
    auto free_params = FreeParameter();
    const auto fn_m_total = TotalDustMassFileName(free_params);
    const auto fn_m       = DustMassDistributionFileName(free_params);
    const auto fn_n       = DustNumberDistributionFileName(free_params);

    asano_model.Calculate(free_params, fn_m_total, fn_m, fn_n);

    return 0;
}

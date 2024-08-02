library(brms)

pop_mod <- brm(seed = 111, save_pars = save_pars(latent=TRUE),
            data=mdat,
            chains=2,control=list(adapt_delta=.995,max_treedepth=15),
            iter=16000,
            Pop.pop.corr. | se(S.E..1, sigma=TRUE) ~ me(xmate_corr, xmate_corr_se) +
            (1|mm(Phenotype.1,Phenotype.2)))

direct_mod <- brm(seed = 111, save_pars = save_pars(latent=TRUE),
            data=mdat,
            chains=2,control=list(adapt_delta=.995,max_treedepth=15),
            iter=16000,
            Direct.direct.corr. | se(S.E., sigma=TRUE) ~ me(xmate_corr, xmate_corr_se) +
            (1|mm(Phenotype.1,Phenotype.2)))

joint_mod <- brm(seed = 111, save_pars = save_pars(latent=TRUE),
            data=mdat,
            chains=2,control=list(adapt_delta=.995,max_treedepth=15),
            iter=16000,
            xmate_corr | se(xmate_corr_se, sigma=TRUE) ~ me(rpop, rpopse)+me(rdir, rdirse) +
            (1|mm(Phenotype.1,Phenotype.2)))

diff_mod <- brm(seed = 111, save_pars = save_pars(latent=TRUE),
            data=mdat,
            chains=2,control=list(adapt_delta=.995,max_treedepth=15),
            iter=16000,
            diff | se(diff_se, sigma=TRUE) ~ me(xmate_corr, xmate_corr_se) +
            (1|mm(Phenotype.1,Phenotype.2)))

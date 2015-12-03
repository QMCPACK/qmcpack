##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_calculations.py                                           #
#    Provides compact interfaces to a few common QMCPACK             #
#    workflows involving PWSCF and SQD.  While this code functions,  #
#    this design route has not been developed further with the       #
#    intent of preserving a common user interface for workflow       #
#    generation regardless of what simulation codes are used.        #
#                                                                    #
#  Content summary:                                                  #
#    basic_qmc                                                       #
#      Function generates a basic QMC workflow including orbital     #
#      generation with PWSCF, orbital conversion with pw2qmcpack,    #
#      Jastrow optimization with QMCPACK, and diffusion Monte Carlo  #
#      with QMCPACK.                                                 #
#                                                                    #
#    standard_qmc                                                    #
#      Same as basic_qmc, except that PWSCF orbital generation       #
#      starts with a charge density calculation (scf) followed by    #
#      non-scf calculations to make the orbitals.                    #
#                                                                    #
#    sqd_qmc                                                         #
#      Function generates a simple QMC workflow for an atom          #
#      starting with Hartree-Fock orbital generation with the SQD    #
#      code and followed by Jastrow optimization and DMC with        #
#      QMCPACK.                                                      #
#                                                                    #
#====================================================================#


import os
from numpy import array,empty

from generic import obj
from machines import Job
from physical_system import generate_physical_system
from pwscf import generate_pwscf
from pw2qmcpack import generate_pw2qmcpack
from pw2casino import generate_pw2casino
from wfconvert import generate_wfconvert
from sqd import generate_sqd
from qmcpack import generate_qmcpack,BundledQmcpack
from qmcpack_input import loop,vmc,dmc,linear,cslinear
from structure import kmesh


def error(msg):
    print 'Error: '+msg
    exit()
#end def error


def sort_pseudos(pseudos):
    dftpseudos = []
    qmcpseudos = []
    for pseudo in pseudos:
        pp = pseudo.lower()
        if pp.endswith('upf') or pp.endswith('ncpp'):
            dftpseudos.append(pseudo)
        elif pp.endswith('xml'):
            qmcpseudos.append(pseudo)
        else:
            error(pseudo+' is not a pwscf or qmcpack pseudopotential')
        #end if
    #end for
    return dftpseudos,qmcpseudos
#end def sort_pseudos


def generate_kpoint(ks,kaxes,grids=True):
    is_string = isinstance(ks,str)
    if is_string:
        kgrid     = 1,1,1
        dft_kgrid = 1,1,1
        if ks=='G':
            kshift    = 0,0,0
            dft_kshift = 0,0,0
        elif ks=='L':
            kshift    = .5,.5,.5
            dft_kshift = 1,1,1
        elif len(ks)>1 and ks[1:].isdigit():
            dft_kshift = tuple(array(tuple(ks[1:]),dtype=int))
            if len(dft_kshift)!=3:
                error('kpoint string must have at least 3 entires for shift\n  you provided {0}\n  input kpoint string: {1}\n  a valid example: {2}'.format(len(dft_kshift),ks,ks[0]+'111'))
            #end if
            kshift = empty((3,),dtype=float)
            for d in range(3):
                dks = dft_kshift[d]
                if dks==1:
                    kshift[d] = .5
                elif dks==0:
                    kshift[d] = 0
                else:
                    error('invalid value in kpoint string\n  valid values are 0 or 1\n  you provided: {0}\n  full input string: {1}'.format(dks,ks))
                #end if
            #end for
            kshift = tuple(kshift)
        #end if
        kpoint = kmesh(kaxes,kgrid,kshift)[0]
    else:
        kpoint = array(ks)
    #end if

    if not grids or not is_string:
        kgrid      = None
        kshift     = None
        dft_kgrid  = None
        dft_kshift = None
    #end if

    k = obj(
        kpoint     = kpoint,
        kgrid      = kgrid,
        kshift     = kshift,
        dft_kgrid  = dft_kgrid,
        dft_kshift = dft_kshift
        )

    return k
#end def generate_kpoint




def basic_qmc(identifier       = '',
              system           = None,
              pseudos          = None,
              functional       = None,
              folded_dft       = True,
              ecut             = 200,
              ecutrho          = None,
              conv_thr         = 1e-8,
              mixing_beta      = .7,
              nosym            = True,
              nscf             = False,
              hubbard_u        = None,
              start_mag        = None,
              kinetic_E        = False, 
              wfconvert        = False,
              bconds           = None,
              jastrows         = 'generateJ12',
              meshfactor       = 1.0,
              det_format       = 'new',
              precision        = 'float',
              vmc_timestep     = None,
              nonlocalmoves    = None,
              perform_opt      = False,
              opt              = None,
              opt_calcs        = None,
              opt_kpoint       = None,
              sample_factor    = 1.0,
              block_opt        = False,
              skip_submit_opt  = False,
              qmc_calcs        = None,
              estimators       = None,
              corrections      = None,
              block_qmc        = False,
              force_write      = False,
              skip_submit_qmc  = False,
              basepath         = '',
              directory        = 'qmc_calc',
              dftdir           = 'scf',
              dftoptdir        = 'scfopt',
              optdir           = 'opt',
              qmcdir           = 'qmc',
              scfjob           = None,
              p2cjob           = None,
              p2qjob           = None,
              wfcjob           = None,
              optjob           = None,
              qmcjob           = None,
              dft_dependencies = None,
              opt_dependencies = None,
              remove_cell      = False,
              return_list      = True):

    if p2cjob is None:
        p2cjob = Job(cores=1)
    #end if
    if p2qjob is None:
        p2qjob = Job(cores=1)
    #end if
    if wfcjob is None:
        wfcjob = Job(cores=1)
    #end if
    if system is None:
        error('system is a required input to basic_qmc')
    #end if
    if pseudos is None:
        error('pseudos is a required input to basic_qmc')
    #end if
    if identifier!='':
        directory = directory+'_'+identifier
        qmcid = identifier
    else:
        qmcid = 'qmc'
    #end if

    dftpath = os.path.join(basepath,directory,dftdir)
    dftoptpath = os.path.join(basepath,directory,dftoptdir)
    optpath = os.path.join(basepath,directory,optdir)
    qmcpath = os.path.join(basepath,directory,qmcdir)

    dftpseudos,qmcpseudos = sort_pseudos(pseudos)

    sys = system.copy()

    structdep = None
    if dft_dependencies!=None:
        for dep in dft_dependencies:
            if dep[1]=='structure':
                structdep = dep
            #end if
        #end for
    #end if

    #add vmc timestep and nonlocalmoves to qmc_calcs, if requested
    if qmc_calcs!=None:
        for calc in qmc_calcs:
            if isinstance(calc,loop):
                calc = calc.qmc
            #end if
            if isinstance(calc,vmc) or isinstance(calc,linear) or isinstance(calc,cslinear):
                if vmc_timestep!=None and not 'timestep' in calc:
                    calc.timestep = vmc_timestep
                #end if
                if nonlocalmoves!=None:
                    calc.nonlocalpp = nonlocalmoves
                #end if
            elif isinstance(calc,dmc):
                if nonlocalmoves!=None:
                    calc.nonlocalmoves = nonlocalmoves
                #end if
            #end if
        #end for
    #end if

                
    if not nscf:
        scf_id = 'scf'
        scf_it = 'scf'
    else:
        scf_id = 'nscf'
        scf_it = 'nscf'
    #end if

    scf = generate_pwscf(
        identifier   = scf_id,
        path         = dftpath,
        job          = scfjob,
        input_type   = scf_it,
        input_dft    = functional,
        ecut         = ecut,
        ecutrho      = ecutrho,
        conv_thr     = conv_thr,
        mixing_beta  = mixing_beta,
        nosym        = nosym,
        hubbard_u    = hubbard_u,
        start_mag    = start_mag,
        pseudos      = dftpseudos,
        system       = sys,
        use_folded   = folded_dft
        )
    p2q = generate_pw2qmcpack(
        identifier   = 'p2q',
        path         = dftpath,
        job          = p2qjob,
        write_psir   = False
        )
    if wfconvert:
        wfc = generate_wfconvert(
            identifier   = 'wfc',
            path         = dftpath,
            job          = wfcjob,
            spline       = False
            )
    #end if
    qmc = generate_qmcpack(
        identifier   = qmcid,
        path         = qmcpath,
        job          = qmcjob,
        block        = block_qmc,
        force_write  = force_write,
        skip_submit  = skip_submit_qmc,
        input_type   = 'basic',
        pseudos      = qmcpseudos,
        system       = sys,
        bconds       = bconds,
        remove_cell  = remove_cell,
        jastrows     = jastrows,
        meshfactor   = meshfactor,
        det_format   = det_format,
        precision    = precision,
        estimators   = estimators,
        corrections  = corrections,
        calculations = qmc_calcs
        )

    sims = obj(
        scf = scf,
        p2q = p2q,
        qmc = qmc
        )

    if dft_dependencies!=None:
        scf.depends(*dft_dependencies)
    #end if
    p2q.depends(scf,'orbitals')
    if wfconvert:
        sims.wfc = wfc
        wfc.depends(p2q,'orbitals')
        qmc.depends(wfc,'orbitals')
    else:
        qmc.depends(p2q,'orbitals')
    #end if

    if structdep!=None:
        qmc.depends(structdep)
    #end if
    if not perform_opt and opt is None and opt_dependencies!=None:
        qmc.depends(*opt_dependencies)
    #end if

    if perform_opt:
        if wfconvert:
            orbdep = wfc
        else:
            orbdep = p2q
        #end if
        opt_sys = sys
        if opt_kpoint!=None:
            kaxes = opt_sys.structure.kaxes.copy()

            k = generate_kpoint(opt_kpoint,kaxes,grids=not nscf)

            opt_kgrid     = k.kgrid
            opt_kshift    = k.kshift
            dftopt_kshift = k.dft_kgrid
            dftopt_kgrid  = k.dft_kshift

            opt_sys = opt_sys.copy()
            opt_sys.structure.clear_kpoints()
            opt_sys.structure.add_kpoints([k.kpoint],[1.0])

            scfopt = generate_pwscf(
                identifier = scf_id,
                path       = dftoptpath,
                job        = scfjob,
                input_type = scf_it,
                input_dft  = functional,
                ecut       = ecut,
                conv_thr   = conv_thr,
                nosym      = nosym,
                hubbard_u  = hubbard_u,
                start_mag  = start_mag,
                pseudos    = dftpseudos,
                system     = opt_sys,
                kgrid      = dftopt_kgrid,
                kshift     = dftopt_kshift
                )
            p2qopt = generate_pw2qmcpack(
                identifier = 'p2q',
                path       = dftoptpath,
                job        = p2qjob,
                write_psir = False
                )
            if wfconvert:
                wfcopt = generate_wfconvert(
                    identifier = 'wfc',
                    path       = dftoptpath,
                    job        = wfcjob,
                    spline     = False
                    )
            #end if
            if dft_dependencies!=None:
                scfopt.depends(*dft_dependencies)
            #end if
            p2qopt.depends(scfopt,'orbitals')
            sims.set(
                scfopt = scfopt,
                p2qopt = p2qopt
                )
            if wfconvert:
                sims.wfcopt = wfcopt
                wfcopt.depends(p2qopt,'orbitals')
                orbdep = wfcopt
            else:
                orbdep = p2qopt
            #end if
        #end if
        optjob.set_processes()
        opt = generate_qmcpack(
            identifier = 'opt',
            path       = optpath,
            job        = optjob,
            block      = block_opt,
            skip_submit= skip_submit_opt,
            force_write= force_write,
            input_type = 'opt_jastrow',
            pseudos    = qmcpseudos,
            system     = opt_sys,
            bconds     = bconds,
            remove_cell= remove_cell,
            jastrows   = jastrows,
            meshfactor = meshfactor,
            det_format = det_format,
            precision  = precision,
            timestep   = vmc_timestep,
            nonlocalpp = nonlocalmoves,
            corrections= [],
            sample_factor = sample_factor,
            opt_calcs  = opt_calcs,
            processes  = optjob.processes,
            threads    = optjob.threads
            )
        opt.depends(orbdep,'orbitals')
        if structdep!=None:
            opt.depends(structdep)
        #end if
        if opt_dependencies!=None:
            opt.depends(*opt_dependencies)
        #end if
        if kinetic_E and opt_kpoint!=None:
            p2copt = generate_pw2casino(
                identifier = 'p2c',
                path       = dftoptpath,
                job        = p2cjob
                )
            p2copt.depends(scfopt,'orbitals')
            sims.p2copt = p2copt
        #end if
    #end if
    if opt is not None:
        sims.opt = opt
        qmc.depends(opt,'jastrow')
    #end if

    if kinetic_E:
        p2c = generate_pw2casino(
            identifier = 'p2c',
            path       = dftpath,
            job        = p2cjob
            )
        p2c.depends(scf,'orbitals')
        sims.p2c = p2c
    #end if

    if return_list:
        order = ['scfopt','p2copt','p2qopt','wfcopt','opt','scf','p2c','p2q','wfc','qmc']
        simlist = []
        for simid in order:
            if simid in sims:
                simlist.append(sims[simid])
                del sims[simid]
            #end if
        #end for
        return simlist
    else:
        return sims
    #end if
#end def basic_qmc




def standard_qmc(identifier      = '',
                 scf_kgrid       = None,
                 scf_kshift      = (1,1,1),
                 system          = None,
                 pseudos         = None,
                 folded_dft      = True,
                 functional      = None,
                 ecut            = 200,
                 ecutrho         = None,
                 conv_thr        = 1e-8,
                 mixing_beta     = .7,
                 nosym           = False,
                 nscf_nosym      = True,
                 hubbard_u       = None,
                 start_mag       = None,
                 kinetic_E       = False, 
                 wfconvert       = False,
                 bconds          = None,
                 jastrows        = 'generateJ12',
                 meshfactor      = 1.0,
                 det_format      = 'new',
                 precision       = 'float',
                 vmc_timestep    = None,
                 nonlocalmoves   = None,
                 perform_opt     = False,  
                 opt             = None,
                 opt_calcs       = None,
                 opt_kpoint      = None,
                 sample_factor   = 1.0,
                 block_opt       = False,
                 skip_submit_opt = False,
                 qmc_calcs       = None,
                 estimators      = None,
                 corrections     = None,
                 block_qmc       = False,
                 force_write     = False,
                 skip_submit_qmc = False,
                 basepath        = '',
                 directory       = 'qmc_calc',
                 scfdir          = 'scf',
                 nscfdir         = 'nscf',
                 nscfoptdir      = 'nscfopt',
                 optdir          = 'opt',
                 qmcdir          = 'qmc',
                 scfjob          = None,
                 nscfjob         = None,
                 p2cjob          = None,
                 p2qjob          = None,
                 wfcjob          = None,
                 optjob          = None,
                 qmcjob          = None,
                 dft_dependencies= None,
                 opt_dependencies= None,
                 remove_cell     = False,
                 return_list     = True):
    if p2cjob is None:
        p2cjob = Job(cores=1)
    #end if
    if p2qjob is None:
        p2qjob = Job(cores=1)
    #end if
    if wfcjob is None:
        wfcjob = Job(cores=1)
    #end if
    if system is None:
        error('system is a required input to standard_qmc')
    #end if
    if pseudos is None:
        error('pseudos is a required input to standard_qmc')
    #end if
    if scf_kgrid is None:
        error('scf_kgrid is a required input to standard_qmc')
    #end if


    scfpath = os.path.join(basepath,directory,scfdir)
    dftpseudos,qmcpseudos = sort_pseudos(pseudos)
    sys = system.copy()    

    scf = generate_pwscf(
        identifier = 'scf',
        path       = scfpath,
        job        = scfjob,
        input_type = 'scf',
        input_dft  = functional,
        ecut       = ecut,
        ecutrho    = ecutrho,
        conv_thr   = conv_thr,
        mixing_beta= mixing_beta,
        nosym      = nosym,
        hubbard_u  = hubbard_u,
        start_mag  = start_mag,
        pseudos    = dftpseudos,
        system     = sys,
        kgrid      = scf_kgrid,
        kshift     = scf_kshift,
        wf_collect = False,
        use_folded = folded_dft
        )
    

    if dft_dependencies is None:
        dft_dependencies = []
    else:
        scf.depends(*dft_dependencies)
        dft_dependencies = list(dft_dependencies)
    #end if
    dft_dependencies.append((scf,'charge_density'))


    sims = basic_qmc(
        nscf             = True            ,   
        identifier       = identifier      ,   
        system           = system          ,   
        pseudos          = pseudos         ,   
        folded_dft       = folded_dft      ,
        functional       = functional      ,   
        jastrows         = jastrows        ,   
        meshfactor       = meshfactor      ,   
        det_format       = det_format      ,
        precision        = precision       ,
        ecut             = ecut            ,   
        ecutrho          = ecutrho         ,   
        conv_thr         = conv_thr        ,   
        mixing_beta      = mixing_beta     ,   
        nosym            = nscf_nosym      ,    
        hubbard_u        = hubbard_u       ,
        start_mag        = start_mag       ,
        kinetic_E        = kinetic_E       ,
        wfconvert        = wfconvert       ,
        bconds           = bconds          ,
        vmc_timestep     = vmc_timestep    ,
        nonlocalmoves    = nonlocalmoves   ,
        perform_opt      = perform_opt     ,   
        opt              = opt             ,   
        opt_calcs        = opt_calcs       ,   
        opt_kpoint       = opt_kpoint      ,   
        sample_factor    = sample_factor   ,   
        block_opt        = block_opt       ,   
        skip_submit_opt  = skip_submit_opt ,
        qmc_calcs        = qmc_calcs       ,   
        estimators       = estimators      ,
        corrections      = corrections     ,
        block_qmc        = block_qmc       ,   
        force_write      = force_write     ,   
        skip_submit_qmc  = skip_submit_qmc ,   
        basepath         = basepath        ,   
        directory        = directory       ,   
        dftdir           = nscfdir         ,   
        dftoptdir        = nscfoptdir      ,   
        optdir           = optdir          ,   
        qmcdir           = qmcdir          ,   
        scfjob           = nscfjob         ,   
        p2cjob           = p2cjob          ,   
        p2qjob           = p2qjob          ,   
        wfcjob           = wfcjob          ,   
        optjob           = optjob          ,   
        qmcjob           = qmcjob          ,   
        dft_dependencies = dft_dependencies,
        opt_dependencies = opt_dependencies,
        remove_cell      = remove_cell,
        return_list      = return_list         
        )

    if return_list:
        sims = [scf]+sims
    else:
        if 'scf' in sims:
            sims.nscf = sims.scf
            del sims.scf
        #end if
        if 'scfopt' in sims:
            sims.nscfopt = sims.scfopt
            del sims.scfopt
        #end if
        sims.scf = scf
    #end if

    return sims
#end def standard_qmc








def sqd_qmc(identifier       = '',
            system           = None,
            filling          = None,
            updown           = None,
            down             = None,
            up               = None,
            grid_type        = 'log',
            ri               = 1e-6,
            rf               = 400,
            npts             = 10001,
            max_iter         = 1000,
            etot_tol         = 1e-8,
            eig_tol          = 1e-12,
            mix_ratio        = 0.7,            
            jastrows         = 'generateJ12',
            sqd_jastrow      = True,
            sqd_jastrow_qmc  = False,
            vmc_timestep     = 0.5,
            perform_opt      = False,
            opt              = None,
            opt_calcs        = None,
            sample_factor    = 1.0,
            block_opt        = False,
            skip_submit_opt  = False,
            qmc_calcs        = None,
            estimators       = None,
            block_qmc        = False,
            force_write      = False,
            skip_submit_qmc  = False,
            basepath         = '',
            directory        = 'qmc_sqd',
            hfdir            = 'hf',
            optdir           = 'opt',
            qmcdir           = 'qmc',
            hfjob            = None,
            optjob           = None,
            qmcjob           = None,
            opt_dependencies = None,
            return_list      = True):
    if system is None:
        error('system is a required input to sqd_qmc')
    #end if
    if identifier!='':
        directory = directory+'_'+identifier
        qmcid = identifier
    else:
        qmcid = 'qmc'
    #end if
    hfid = 'hf'

    hfpath = os.path.join(basepath,directory,hfdir)
    optpath = os.path.join(basepath,directory,optdir)
    qmcpath = os.path.join(basepath,directory,qmcdir)

    sys = system.copy()


    #add vmc timestep and nonlocalmoves to qmc_calcs, if requested
    if qmc_calcs!=None and vmc_timestep!=None:
        for calc in qmc_calcs:
            if isinstance(calc,loop):
                calc = calc.qmc
            #end if
            if isinstance(calc,vmc) or isinstance(calc,linear) or isinstance(calc,cslinear) and not 'timestep' in calc:
                calc.timestep = vmc_timestep
            #end if
        #end for
    #end if



    hf = generate_sqd(
        identifier = hfid,
        path       = hfpath,
        job        = hfjob,
        system     = sys,
        filling    = filling,
        updown     = updown,
        up         = up,
        down       = down,
        grid_type  = grid_type  ,
        ri         = ri         ,
        rf         = rf         ,
        npts       = npts       ,
        max_iter   = max_iter   ,
        etot_tol   = etot_tol   ,
        eig_tol    = eig_tol    ,
        mix_ratio  = mix_ratio  
        )
    qmc = generate_qmcpack(
        identifier = qmcid,
        path       = qmcpath,
        job        = qmcjob,
        block      = block_qmc,
        force_write= force_write,
        skip_submit= skip_submit_qmc,
        input_type = 'basic',
        system     = sys,
        jastrows   = jastrows,
        estimators = estimators,
        calculations = qmc_calcs
        )
    qmc.depends(hf,'orbitals')
    if sqd_jastrow_qmc:
        qmc.depends(hf,'jastrow')
    #end if

    sims = obj(
        hf = hf,
        qmc = qmc
        )

    if perform_opt:
        opt_sys = sys
        optjob.set_processes()
        opt = generate_qmcpack(
            identifier = 'opt',
            path       = optpath,
            job        = optjob,
            block      = block_opt,
            skip_submit= skip_submit_opt,
            force_write= force_write,
            input_type = 'opt_jastrow',
            system     = opt_sys,
            jastrows   = jastrows,
            timestep   = vmc_timestep,
            sample_factor = sample_factor,
            opt_calcs  = opt_calcs,
            processes  = optjob.processes,
            threads    = optjob.threads
            )
    #end if
    if opt is not None:
        if sqd_jastrow:
            opt.depends(hf,'jastrow')
        #end if
        opt.depends(hf,'orbitals')
        if opt_dependencies!=None:
            opt.depends(*opt_dependencies)
        #end if
        sims.opt = opt
        qmc.depends(opt,'jastrow')
    #end if

    if return_list:
        order = ['hf','opt','qmc']
        simlist = []
        for simid in order:
            if simid in sims:
                simlist.append(sims[simid])
                del sims[simid]
            #end if
        #end for
        return simlist
    else:
        return sims
    #end if
#end def sqd_qmc










################################################################
#                                                              #
#      These are antiquated, do not use in present form!       #
#                                                              #
################################################################


def optimization(identifier = '',
                 system     = None,
                 pseudos    = None,
                 functional = None,
                 jastrows   = 'generateJ12',
                 meshfactor = 1.0,
                 ecut       = 200,
                 conv_thr   = 1e-8,
                 kinetic_E  = True,
                 block_qmc  = False,
                 force_write= False,
                 opt_calcs  = None,
                 opt_kpoint = None,
                 basepath   = '',
                 directory  = 'optimization',
                 dftdir     = 'dft',
                 optdir     = 'opt',
                 dftjob     = None,
                 p2cjob     = None,
                 p2qjob     = None,
                 wfcjob     = None,
                 optjob     = None):
    error('do not use this function (optimization) unless you update it!')
    if p2cjob is None:
        p2cjob = Job(cores=1)
    #end if
    if p2qjob is None:
        p2qjob = Job(cores=1)
    #end if
    if wfcjob is None:
        wfcjob = Job(cores=1)
    #end if
    if system is None:
        error('system is a required input to optimization')
    #end if
    if pseudos is None:
        error('pseudos is a required input to optimization')
    #end if
    if identifier!='':
        directory = directory+'_'+identifier
        qmcid = identifier
    else:
        qmcid = 'qmc'
    #end if

    dftpath = os.path.join(basepath,directory,dftdir)
    optpath = os.path.join(basepath,directory,optdir)

    dftpseudos,qmcpseudos = sort_pseudos(pseudos)

    sys = system.copy()
    if opt_kpoint!=None:
        sys.structure.kpoints = array([array(opt_kpoint)])
    #end if
    optjob.set_processes()

    if sys.folded_system!=None:
        dftsys = sys.folded_system
    else:
        dftsys = sys
    #end if

    scf = generate_pwscf(
        identifier = 'scf',
        path       = dftpath,
        job        = dftjob,
        input_type = 'scf',
        input_dft  = functional,
        ecut       = ecut,
        conv_thr   = conv_thr,
        pseudos    = dftpseudos,
        system     = dftsys
        )
    p2q = generate_pw2qmcpack(
        identifier = 'p2q',
        path       = dftpath,
        job        = p2qjob,
        write_psir = False
        )
    wfc = generate_wfconvert(
        identifier = 'wfc',
        path       = dftpath,
        job        = wfcjob,
        spline     = False
        )
    opt = generate_qmcpack(
        identifier = 'opt',
        path       = optpath,
        job        = optjob,
        block      = block_qmc,
        force_write= force_write,
        input_type = 'opt_jastrow',
        pseudos    = qmcpseudos,
        system     = sys,
        jastrows   = jastrows,
        meshfactor = meshfactor,
        opt_calcs  = opt_calcs,
        processes  = optjob.processes,
        threads    = optjob.threads
        )
    p2q.depends(scf,'orbitals')
    wfc.depends(p2q,'orbitals')
    opt.depends(wfc,'orbitals')

    if kinetic_E:
        p2c = generate_pw2casino(
            identifier = 'p2c',
            path       = dftpath,
            job        = p2cjob
            )
        p2c.depends(scf,'orbitals')
        sims = [scf,p2c,p2q,wfc,opt]        
    else:
        sims = [scf,p2q,wfc,opt]
    #end if

    return sims
#end def optimization

    


def twist_averaged_qmc(identifier = '',
                       system     = None,
                       pseudos    = None,
                       functional = None,
                       jastrows   = 'generateJ12',
                       meshfactor = 1.0,
                       ecut       = 200,
                       conv_thr   = 1e-8,
                       kgrid      = None,
                       kshift     = None,
                       kinetic_E  = True,
                       shared_orb = True,
                       perform_opt= False,
                       opt        = None,
                       opt_calcs  = None,
                       opt_kpoint = (0,0,0),
                       block_opt  = False,
                       qmc_calcs  = None,
                       force_write= False,
                       block_qmc  = False,
                       basepath   = '',
                       directory  = 'twist_average',
                       dftdir     = 'dft',
                       optdir     = 'opt',
                       qmcdir     = 'qmc',
                       dftjob     = None,
                       p2cjob     = None,
                       p2qjob     = None,
                       wfcjob     = None,
                       optjob     = None,
                       qmcjob     = None):
    error('do not use this function (twist_averaged_qmc) unless you update it!')
    if p2cjob is None:
        p2cjob = Job(cores=1)
    #end if
    if p2qjob is None:
        p2qjob = Job(cores=1)
    #end if
    if wfcjob is None:
        wfcjob = Job(cores=1)
    #end if

    sims = []
    
    if system is None:
        error('system is a required input to basic_qmc')
    #end if
    if pseudos is None:
        error('pseudos is a required input to basic_qmc')
    #end if
    if identifier!='':
        directory = directory+'_'+identifier
        qmcid = identifier
    else:
        qmcid = 'qmc'
    #end if

    dftpseudos,qmcpseudos = sort_pseudos(pseudos)

    structure = system.structure
    if kgrid!=None and kshift!=None:
        system = system.copy()
        structure = system.structure
        structure.clear_kpoints()
        structure.add_kmesh(kgrid,kshift)
    #end if
    if shared_orb:
        kpoints = ['shared_orb']
    else:
        kpoints = structure.kpoints
    #end if

    optpath = os.path.join(basepath,directory,optdir)
    qmcpath = os.path.join(basepath,directory,qmcdir)

    if perform_opt and opt is None:
        optsims = optimization(
            identifier = identifier,
            system     = system,
            pseudos    = pseudos,
            functional = functional,
            jastrows   = jastrows,
            meshfactor = meshfactor,
            ecut       = ecut,
            conv_thr   = conv_thr,
            kinetic_E  = kinetic_E,
            block_qmc  = block_opt,
            force_write= force_write,
            opt_calcs  = opt_calcs,
            opt_kpoint = opt_kpoint,
            basepath   = os.path.join(basepath,directory),
            dftjob     = dftjob,
            p2cjob     = p2cjob,
            p2qjob     = p2qjob,
            wfcjob     = wfcjob,
            optjob     = optjob       
            )
        opt = optsims[-1]
        sims.extend(optsims)
    #end if

    n=0
    qmcs = []
    for kpoint in kpoints:
        if shared_orb:
            sys = system
            idext = ''
        else:
            sys = system.copy()
            sys.structure.kpoints = array([kpoint])
            idext = '_'+str(n).zfill(3)
        #end if

        dftpath = os.path.join(basepath,directory,dftdir+idext)

        if sys.folded_system!=None:
            dftsys = sys.folded_system
        else:
            dftsys = sys
        #end if


        scf = generate_pwscf(
            identifier = 'scf',
            path       = dftpath,
            job        = dftjob,
            input_type = 'scf',
            input_dft  = functional,
            ecut       = ecut,
            conv_thr   = conv_thr,
            pseudos    = dftpseudos,
            system     = dftsys
            )
        p2q = generate_pw2qmcpack(
            identifier = 'p2q',
            path       = dftpath,
            job        = p2qjob,
            write_psir = False
            )
        wfc = generate_wfconvert(
            identifier = 'wfc',
            path       = dftpath,
            job        = wfcjob,
            spline     = False
            )
        qmc = generate_qmcpack(
            identifier = qmcid+'_twist'+idext,
            path       = qmcpath,
            job        = qmcjob,
            block      = block_qmc,
            force_write= force_write,
            input_type = 'basic',
            pseudos    = qmcpseudos,
            system     = sys,
            jastrows   = jastrows,
            meshfactor = meshfactor,
            calculations = qmc_calcs
            )
        
        qmcs.append(qmc)

        p2q.depends(scf,'orbitals')
        wfc.depends(p2q,'orbitals')
        qmc.depends(wfc,'orbitals')
        if opt is not None:
            qmc.depends(opt,'jastrow')
        #end if
        sims += [scf,p2q,wfc]

        if kinetic_E:
            p2c = generate_pw2casino(
                identifier = 'p2c',
                path       = dftpath,
                job        = p2cjob
                )
            p2c.depends(scf,'orbitals')
            sims.append(p2c)
        #end if

        n+=1
    #end for

    if shared_orb:
        sims.extend(qmcs)
    else:
        taqmc = BundledQmcpack(
            identifier = identifier,
            path       = qmcpath,
            job        = qmcjob,
            block      = block_qmc,
            force_write= force_write,
            app_name   = 'qmcapp_complex',
            sims       = qmcs,
            system     = system
        )        
        sims.append(taqmc)
        for qmc in qmcs:
            qmc.block = True
        #end for
    #end if

    return sims
#end def twist_averaged_qmc
    



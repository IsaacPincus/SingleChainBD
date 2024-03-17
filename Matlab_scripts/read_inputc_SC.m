function simulation_parameter_structure = read_inputc_SC()
    Stype = 1;
    NBeads = 5;
    hstar = 0;
    EV = 0;
    zstar = 0;
    dstar = 1;
    MaxGamma = 0;
    LookupTableOpt = 0;
    LookupTableTolerance = 0;
    BendingPotentialType = 0;
    BendingStiffness = 0;
    NaturalAngleData = 0;
    NaturalAngleFromFile = 0;
    sqrtb = 7;
    Q0s = 1;
    gdots = 0;
    emax = 0;
    bens = 0;
    Netcdf = 0;
    Restart = 0;
    Nsamples = 10;
    phiFromFile = 0;
    phi_data = 0;
    Teq = 3;
    Tpr = 5;
    ntrajinput = 0;
    ndelts = 1;
    variance_reduction = 0;
    FlowTypeProduction = 1;
    InitialConfiguration = 1;
    COM_update_on = 1;
    TrapOneStrength = 0;
    TrapTwoStrength = 0;
    TrapOneInitialPosition = 0;
    TrapTwoInitialPosition = 0;
    TrapTwoVelocity = 0;
        

    inputc_table = readtable('inputc.dat');
    param_names = inputc_table.Var1;
    param_values = inputc_table.Var2;

    count = 0;
    for i = 1:length(param_names)
        param_name = param_names{i};
        param_value = param_values(i);
        switch (param_name)
            case {"SpType","SType"}
                Stype = param_value;
            case {"COM_update_on"}
                COM_update_on = param_value;
            case ("NBeads")
                NBeads = param_value;
            case ("hstar")
                hstar = param_value;
            case ("EV")
                EV = param_value;
            case {"zstarByEpsstar","zstar"}
                zstar = param_value;
            case ("dstar")
                dstar = param_value;
            case ("MaxGamma")
                MaxGamma = param_value;
            case ("LookupTableOpt")
                LookupTableOpt = param_value;
            case ("LookupTableTolerance")
                LookupTableTolerance = param_value;
            case {"BendingPotential", "bending_potential",...
                    "BendingPotentialType"}
                BendingPotentialType = param_value;
            case {"BendingStiffness", "bending_stiffness"}
                BendingStiffness = param_value;
            case {"natural_angle_data", "natural_angle_scalar_input"}
                NaturalAngleData = param_value;
            case {"natural_angle_from_file", "NaturalAngleFromFile"}
                NaturalAngleFromFile = param_value;
            case {"TrapOneInitialPosition"}
                TrapOneInitialPosition = param_value;
            case {"TrapTwoInitialPosition"}
                TrapTwoInitialPosition = param_value;
            case {"TrapOneStrength"}
                TrapOneStrength = param_value;
            case {"TrapTwoStrength"}
                TrapTwoStrength = param_value;
            case {"TrapTwoVelocity"}
                TrapTwoVelocity = param_value;
            case {"L0star", "dQ", "sqrtb"}
                sqrtb = param_value;
            case {"Q0star","Q0s","sigma"}
                Q0s = param_value;
            case {"Gdot","gdots"}
                gdots = param_value;
            case {"Bens","bens"}
                bens = param_value;
            case {"NetCDF","Netcd"}
                Netcdf = param_value;
            case ("Restart")
                Restart = param_value;
            case ("emax")
                emax = param_value;
            case {"Nsamp","Nsamples"}
                Nsamples = param_value;
            case {"phiFile","phiFromFile"}
                phiFromFile = param_value;
            case {"phiSDK","phi_data"}
                phi_data = param_value;
            case ("Teq")
                Teq = param_value;
            case ("Tpr")
                Tpr = param_value;
            case {"Ntrajdone","ntrajinput"}
                ntrajinput = param_value;
            case ("ndelts")
                ndelts = param_value;
            case ("dtseq")
                count=count+1;
                dtseq(count) = param_value;
            case ("dtsne")
                dtsne(count) = param_value;
            case ("nblock")
                nblock(count) = param_value;
            case ("ntot")
                ntot(count) = param_value;
            case ("tol")
                tol(count) = param_value;
            case ("variance_reduction")
                variance_reduction = param_value;
            case ("FlowType")
                FlowTypeProduction = param_value;
            case {"InitialConfiguration", "Initial_config"}
                InitialConfiguration = param_value;
            otherwise
                disp(['variable not recognised, var= ', param_name]);
        end
    end

    % Create a structure containing the relevant simulation parameters above
    simulation_parameter_structure = struct('Stype', Stype);
    simulation_parameter_structure.COM_update_on = COM_update_on;
    simulation_parameter_structure.NBeads = NBeads;
    simulation_parameter_structure.hstar = hstar;
    simulation_parameter_structure.EV = EV;
    simulation_parameter_structure.zstar = zstar;
    simulation_parameter_structure.dstar = dstar;
    simulation_parameter_structure.MaxGamma = MaxGamma;
    simulation_parameter_structure.LookupTableOpt = LookupTableOpt;
    simulation_parameter_structure.LookupTableTolerance = LookupTableTolerance;
    simulation_parameter_structure.BendingPotentialType = BendingPotentialType;
    simulation_parameter_structure.BendingStiffness = BendingStiffness;
    simulation_parameter_structure.NaturalAngleData = NaturalAngleData;
    simulation_parameter_structure.NaturalAngleFromFile = NaturalAngleFromFile;
    simulation_parameter_structure.TrapOneStrength = TrapOneStrength;
    simulation_parameter_structure.TrapTwoStrength = TrapTwoStrength;
    simulation_parameter_structure.TrapOneInitialPosition = TrapOneInitialPosition;
    simulation_parameter_structure.TrapTwoInitialPosition = TrapTwoInitialPosition;
    simulation_parameter_structure.TrapTwoVelocity = TrapTwoVelocity;
    simulation_parameter_structure.sqrtb = sqrtb;
    simulation_parameter_structure.Q0s = Q0s;
    simulation_parameter_structure.gdots = gdots;
    simulation_parameter_structure.emax = emax;
    simulation_parameter_structure.bens = bens;
    simulation_parameter_structure.Netcdf = Netcdf;
    simulation_parameter_structure.Nsamples = Nsamples;
    simulation_parameter_structure.phiFromFile = phiFromFile;
    simulation_parameter_structure.phi_data = phi_data;
    simulation_parameter_structure.Teq = Teq;
    simulation_parameter_structure.Tpr = Tpr;
    simulation_parameter_structure.ntrajinput = ntrajinput;
    simulation_parameter_structure.ndelts = ndelts;
    simulation_parameter_structure.variance_reduction = variance_reduction;
    simulation_parameter_structure.FlowTypeProduction = FlowTypeProduction;
    simulation_parameter_structure.dtseq = dtseq;
    simulation_parameter_structure.dtsne = dtsne;
    simulation_parameter_structure.nblock = nblock;
    simulation_parameter_structure.ntot = ntot;
    simulation_parameter_structure.tol = tol;
    simulation_parameter_structure.InitialConfiguration = InitialConfiguration;
    simulation_parameter_structure.Restart = Restart;

end
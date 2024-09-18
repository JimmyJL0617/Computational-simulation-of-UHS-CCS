# -*- coding: utf-8 -*-
import os

openfoam_header = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://www.openfoam.com
    \\  /    A nd           | Version:  9
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {file_class};  
    location    "{location}";  
    object      {file_object}; 
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

file_path = 'E:\\generated_openfoampar\\suizhihui\\'

if not os.path.exists(file_path):
    os.makedirs(file_path)
'''
def CO2():
    Total_patch = 132
    Inlet_patch = 72
    Outlet_patch = 104
    up = None
    down = None
    CA_CO2 = 100
    CA_water = 80
    compressible = False

    file = file_path + "patches.txt"
    with open(file, 'w') as f:
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name Inlet;\n')
                f.write('        }\n')
            elif i == Outlet_patch:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name Outlet;\n')
                f.write('        }\n')
            elif i == up:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name up;\n')
                f.write('        }\n')
            elif i == down:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name down;\n')
                f.write('        }\n')
            else:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name wall{};\n'.format(i))
                f.write('        }\n')

    file_U = file_path + "U.txt"
    with open(file_U, 'w') as f:
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            codedFixedValue;\n')
                #f.write('        inletvalue            uniform    (0 0 0);\n')
                f.write('        value            uniform    (0 0 0);\n')
                f.write('        name             inletVelocity;\n')
                f.write('        code\n')
                f.write('        #{\n')
                f.write('            const double t = this->db().time().value();\n')
                f.write('            operator == (vector(0.02*((1/(1+exp(-1e6*t)))-0.5), 0, 0));\n')
                f.write('         #};\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            inletOutlet;\n')
                f.write('        inletValue      uniform (0 0 0);\n')
                f.write('        value           uniform (0 0 0);\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            noSlip;\n')
                f.write('    }\n')
                
    file_p = file_path + "p.txt"
    with open(file_p, 'w') as f:
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            calculated;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            calculated;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            calculated;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')

    file_p_rgh = file_path + "p_rgh.txt"
    with open(file_p_rgh, 'w') as f:
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                #f.write('        inletvalue            uniform    1;\n')
                #f.write('        value            uniform    10.005e6;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value           uniform 0;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                f.write('    }\n')            

    file_T = file_path + "T.txt"
    with open(file_T, 'w') as f:
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    308;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n') 
                
    file_alpha_CO2 =  file_path + "alpha_CO2.txt"  
    with open(file_alpha_CO2, 'w') as f:
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    1;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            zeroGradient;\n')
                #f.write('        thetaProperties            ((CO2 water)  135 0 0 0);\n')
                #f.write('        limit            gradient;\n')
                #f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                if compressible:
                    f.write('        type            alphaContactAngle;\n')
                    #f.write('        theta0            {};\n'.format(CA_CO2))
                    f.write('        thetaProperties   ((H2 water) {} 0 0 0);\n'.format(CA_CO2))
                    f.write('        value            uniform    0;\n')
                    f.write('    }\n')
                else:
                    f.write('        type            constantAlphaContactAngle;\n')
                    f.write('        theta0            {};\n'.format(CA_CO2))
                    f.write('        limit            gradient;\n')
                    f.write('        value            uniform    0;\n')
                    f.write('    }\n')
            
    file_alpha_water =  file_path + "alpha_water.txt"  
    with open(file_alpha_water, 'w') as f:
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            zeroGradient;\n')
                #f.write('        thetaProperties            ((CO2 water)  135 0 0 0);\n')
                #f.write('        limit            gradient;\n')
                #f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                if compressible:
                    f.write('        type            alphaContactAngle;\n')
                    #f.write('        theta0            {};\n'.format(CA_CO2))
                    f.write('        thetaProperties   ((water H2) {} 0 0 0);\n'.format(CA_water))
                    f.write('        value            uniform    0;\n')
                    f.write('    }\n')
                else:
                    f.write('        type            constantAlphaContactAngle;\n')
                    f.write('        theta0            {};\n'.format(CA_water))
                    f.write('        limit            gradient;\n')
                    f.write('        value            uniform    0;\n')
                    f.write('    }\n')
'''
    
def CO2(file_path, Total_patch, Inlet_patch, Outlet_patch, up, down, front, back, CA_CO2, CA_water, compressible = False):
    if not os.path.exists(file_path):
        os.makedirs(file_path)   

    file = file_path + "patches.txt"
    with open(file, 'w') as f:
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name Inlet;\n')
                f.write('        }\n')
            elif i == Outlet_patch:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name Outlet;\n')
                f.write('        }\n')
            elif i == up:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name up;\n')
                f.write('        }\n')
            elif i == down:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name down;\n')
                f.write('        }\n')
            elif i == front:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name front;\n')
                f.write('        }\n')
            elif i == back:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name back;\n')
                f.write('        }\n')
            else:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name wall{};\n'.format(i))
                f.write('        }\n')

    file_U = file_path + "U"
    header = openfoam_header.format(file_class="volVectorField", location="0", file_object="U")
    with open(file_U, 'w') as f:
        f.write(header)
        f.write('dimensions    [0 1 -1 0 0 0 0];\n\n')
        f.write('internalField    uniform (0 0 0);\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            codedFixedValue;\n')
                #f.write('        inletvalue            uniform    (0 0 0);\n')
                f.write('        value            uniform    (0 0 0);\n')
                f.write('        name             inletVelocity;\n')
                f.write('        code\n')
                f.write('        #{\n')
                f.write('            const double t = this->db().time().value();\n')
                f.write('            operator == (vector(0.02*((1/(1+exp(-1e6*t)))-0.5), 0, 0));\n')
                f.write('         #};\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            inletOutlet;\n')
                f.write('        inletValue      uniform (0 0 0);\n')
                f.write('        value           uniform (0 0 0);\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            noSlip;\n')
                f.write('    }\n')
        f.write('}')
                
    file_p = file_path + "p"
    header = openfoam_header.format(file_class="volScalarField", location="0", file_object="p")
    with open(file_p, 'w') as f:
        f.write(header)
        f.write('dimensions    [1 -1 -2 0 0 0 0];\n\n')
        f.write('internalField    uniform 0;\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            calculated;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            calculated;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            calculated;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
        f.write('}')

    file_p_rgh = file_path + "p_rgh"
    header = openfoam_header.format(file_class="volScalarField", location="0", file_object="p_rgh")
    with open(file_p_rgh, 'w') as f:
        f.write(header)
        f.write('dimensions    [1 -1 -2 0 0 0 0];\n\n')
        f.write('internalField    uniform 0;\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                #f.write('        inletvalue            uniform    1;\n')
                #f.write('        value            uniform    10.005e6;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value           uniform 1e7;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                f.write('    }\n')   
        f.write('}')
                
    file_alpha_CO2 =  file_path + "alpha.CO2.orig" 
    header = openfoam_header.format(file_class="volScalarField", location="0", file_object="alpha_CO2")
    with open(file_alpha_CO2, 'w') as f:
        f.write(header)
        f.write('dimensions    [0 0 0 0 0 0 0];\n\n')
        f.write('internalField    uniform 0;\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    1;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            zeroGradient;\n')
                #f.write('        thetaProperties            ((CO2 water)  135 0 0 0);\n')
                #f.write('        limit            gradient;\n')
                #f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            constantAlphaContactAngle;\n')
                f.write('        theta0            {};\n'.format(CA_CO2))
                f.write('        limit            gradient;\n')
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
        f.write('}')
            
    file_alpha_water =  file_path + "alpha.water.orig" 
    header = openfoam_header.format(file_class="volScalarField", location="0", file_object="alpha_water")
    with open(file_alpha_water, 'w') as f:
        f.write(header)
        f.write('dimensions    [0 0 0 0 0 0 0];\n\n')
        f.write('internalField    uniform 1;\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            zeroGradient;\n')
                #f.write('        thetaProperties            ((CO2 water)  135 0 0 0);\n')
                #f.write('        limit            gradient;\n')
                #f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            constantAlphaContactAngle;\n')
                f.write('        theta0            {};\n'.format(CA_CO2))
                f.write('        limit            gradient;\n')
                #f.write('        thetaProperties   ((water H2) {} 0 0 0);\n'.format(CA_water))
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
        f.write('}')
                
def VE():
    Total_patch = 15
    HPAM = 1
    WO = 4
    NO = 10
    outlet = 7
    up = None
    down = None
    CA_HPAM = 60
    CA_oil = 120
    compressible = False

    file = file_path + "patches.txt"
    with open(file, 'w') as f:
        for i in range(Total_patch):
            if i == HPAM:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name HPAM;\n')
                f.write('        }\n')
            elif i == outlet:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name outlet;\n')
                f.write('        }\n')
            elif i == WO:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name WO;\n')
                f.write('        }\n')
            elif i == NO:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name NO;\n')
                f.write('        }\n')
            else:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name wall{};\n'.format(i))
                f.write('        }\n')

    file_U = file_path + "U.txt"
    with open(file_U, 'w') as f:
        for i in range(Total_patch):
            if i == HPAM:
                f.write('    HPAM\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        inletvalue            uniform    (0.02 0 0);\n')
                #f.write('        value            uniform    (0 0 0);\n')
                #f.write('        name             inletVelocity;\n')
                #f.write('        code\n')
                #f.write('        #{\n')
                #f.write('            const double t = this->db().time().value();\n')
                #f.write('            operator == (vector(0.02*((1/(1+exp(-1e6*t)))-0.5), 0, 0));\n')
                #f.write('         #};\n')
                f.write('    }\n')
            elif i == outlet:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            inletOutlet;\n')
                f.write('        inletValue      uniform (0 0 0);\n')
                f.write('        value           uniform (0 0 0);\n')
                f.write('    }\n')
            elif i == WO:
                f.write('    WO\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        inletvalue            uniform    (0 -0.015 0);\n')
                f.write('    }\n')
            elif i == NO:
                f.write('    NO\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        inletvalue            uniform    (0 0.015 0);\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            noSlip;\n')
                f.write('    }\n')
                
    file_p = file_path + "p.txt"
    with open(file_p, 'w') as f:
        for i in range(Total_patch):
            if i == HPAM:
                f.write('    HPAM\n')
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                #f.write('        value           $internalField;\n')
                f.write('    }\n')
            elif i == outlet:
                f.write('    outlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value           uniform  0;\n')
                f.write('    }\n')
            elif i == WO:
                f.write('    WO\n')
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                f.write('    }\n')
            elif i == NO:
                f.write('    NO\n')
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                #f.write('        value           $internalField;\n')
                f.write('    }\n')

                
    file_alpha_HPAM =  file_path + "alpha_HPAM.txt"  
    with open(file_alpha_HPAM, 'w') as f:
        for i in range(Total_patch):
            if i == HPAM:
                f.write('    HPAM\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    1;\n')
                f.write('    }\n')
            elif i == outlet:
                f.write('    outlet\n')
                f.write('    {\n')
                f.write('        type            zeroGradient;\n')
                #f.write('        thetaProperties            ((CO2 water)  135 0 0 0);\n')
                #f.write('        limit            gradient;\n')
                #f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == WO:
                f.write('    WO\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == NO:
                f.write('    NO\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                if compressible:
                    f.write('        type            alphaContactAngle;\n')
                    #f.write('        theta0            {};\n'.format(CA_CO2))
                    f.write('        thetaProperties   ((HPAM oil) {} 0 0 0);\n'.format(CA_HPAM))
                    f.write('        value            uniform    0;\n')
                    f.write('    }\n')
                else:
                    f.write('        type            constantAlphaContactAngle;\n')
                    f.write('        theta0            {};\n'.format(CA_HPAM))
                    f.write('        limit            gradient;\n')
                    f.write('        value            uniform    0;\n')
                    f.write('    }\n')
            
    file_alpha_oil =  file_path + "alpha_oil.txt"  
    with open(file_alpha_oil, 'w') as f:
        for i in range(Total_patch):
            if i == HPAM:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == outlet:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            zeroGradient;\n')
                #f.write('        thetaProperties            ((CO2 water)  135 0 0 0);\n')
                #f.write('        limit            gradient;\n')
                #f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == WO:
                f.write('    WO\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value            uniform    1;\n')
                f.write('    }\n')
            elif i == NO:
                f.write('    NO\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value            uniform    1;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                if compressible:
                    f.write('        type            alphaContactAngle;\n')
                    #f.write('        theta0            {};\n'.format(CA_CO2))
                    f.write('        thetaProperties   ((oil HPAM) {} 0 0 0);\n'.format(CA_oil))
                    f.write('        value            uniform    0;\n')
                    f.write('    }\n')
                else:
                    f.write('        type            constantAlphaContactAngle;\n')
                    f.write('        theta0            {};\n'.format(CA_oil))
                    f.write('        limit            gradient;\n')
                    f.write('        value            uniform    0;\n')
                    f.write('    }\n')
                    
def H2(file_path, Total_patch, Inlet_patch, Outlet_patch, up, down, CA_H2, CA_water, compressible = True):
    if not os.path.exists(file_path):
        os.makedirs(file_path)   

    file = file_path + "patches.txt"
    with open(file, 'w') as f:
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name Inlet;\n')
                f.write('        }\n')
            elif i == Outlet_patch:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name Outlet;\n')
                f.write('        }\n')
            elif i == up:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name up;\n')
                f.write('        }\n')
            elif i == down:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name down;\n')
                f.write('        }\n')
            else:
                f.write('        patch{}\n'.format(i))
                f.write('        {\n')
                f.write('            name wall{};\n'.format(i))
                f.write('        }\n')

    file_U = file_path + "U"
    header = openfoam_header.format(file_class="volVectorField", location="0", file_object="U")
    with open(file_U, 'w') as f:
        f.write(header)
        f.write('dimensions    [0 1 -1 0 0 0 0];\n\n')
        f.write('internalField    uniform (0 0 0);\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                #f.write('        type            codedFixedValue;\n')
                #f.write('        inletvalue            uniform    (0 0 0);\n')
                f.write('        type             fixedValue;\n')
                f.write('        value            uniform    (0.01 0 0);\n')
                #f.write('        name             inletVelocity;\n')
                #f.write('        code\n')
                #f.write('        #{\n')
                #f.write('            const double t = this->db().time().value();\n')
                #f.write('            operator == (vector(0.004*((1/(1+exp(-1e6*t)))-0.5), 0, 0));\n')
                #f.write('         #};\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            pressureInletOutletVelocity;\n')
                #f.write('        inletValue      uniform (0 0 0);\n')
                f.write('        value           uniform (0 0 0);\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            noSlip;\n')
                f.write('    }\n')
        f.write('}')
                
    file_p = file_path + "p"
    header = openfoam_header.format(file_class="volScalarField", location="0", file_object="p")
    with open(file_p, 'w') as f:
        f.write(header)
        f.write('dimensions    [1 -1 -2 0 0 0 0];\n\n')
        f.write('internalField    uniform 1.5e7;\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            calculated;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            calculated;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            calculated;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
        f.write('}')

    file_p_rgh = file_path + "p_rgh"
    header = openfoam_header.format(file_class="volScalarField", location="0", file_object="p_rgh")
    with open(file_p_rgh, 'w') as f:
        f.write(header)
        f.write('dimensions    [1 -1 -2 0 0 0 0];\n\n')
        f.write('internalField    uniform 1.5e7;\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    1.5e7;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            totalPressure;\n')
                f.write('        p0              uniform 1.5e7;\n')
                f.write('        rho             thermo:rho;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            fixedFluxPressure;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')   
        f.write('}')

    file_T = file_path + "T"
    header = openfoam_header.format(file_class="volScalarField", location="0", file_object="T")
    with open(file_T, 'w') as f:
        f.write(header)
        f.write('dimensions    [0 0 0 1 0 0 0];\n\n')
        f.write('internalField    uniform 318;\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    318;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                f.write('        value           $internalField;\n')
                f.write('    }\n')
        f.write('}')
                
    file_alpha_H2 =  file_path + "alpha.H2.orig" 
    header = openfoam_header.format(file_class="volScalarField", location="0", file_object="alpha_H2")
    with open(file_alpha_H2, 'w') as f:
        f.write(header)
        f.write('dimensions    [0 0 0 0 0 0 0];\n\n')
        f.write('internalField    uniform 0;\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    1;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            zeroGradient;\n')
                #f.write('        thetaProperties            ((CO2 water)  135 0 0 0);\n')
                #f.write('        limit            gradient;\n')
                #f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            constantAlphaContactAngle;\n')
                f.write('        theta0            {};\n'.format(CA_H2))
                f.write('        limit            gradient;\n')
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
        f.write('}')
            
    file_alpha_water =  file_path + "alpha.water.orig" 
    header = openfoam_header.format(file_class="volScalarField", location="0", file_object="alpha_water")
    with open(file_alpha_water, 'w') as f:
        f.write(header)
        f.write('dimensions    [0 0 0 0 0 0 0];\n\n')
        f.write('internalField    uniform 1;\n\n')
        f.write('boundaryField\n\n')
        f.write('{\n')
        for i in range(Total_patch):
            if i == Inlet_patch:
                f.write('    Inlet\n')
                f.write('    {\n')
                f.write('        type            fixedValue;\n')
                #f.write('        inletvalue            uniform    1;\n')
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == Outlet_patch:
                f.write('    Outlet\n')
                f.write('    {\n')
                f.write('        type            zeroGradient;\n')
                #f.write('        thetaProperties            ((CO2 water)  135 0 0 0);\n')
                #f.write('        limit            gradient;\n')
                #f.write('        value            uniform    0;\n')
                f.write('    }\n')
            elif i == up:
                f.write('    up\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            elif i == down:
                f.write('    down\n')
                f.write('    {\n')
                f.write('        type            empty;\n')
                f.write('    }\n')
            else:
                f.write('    wall{}\n'.format(i))
                f.write('    {\n')
                f.write('        type            constantAlphaContactAngle;\n')
                #f.write('        theta0            {};\n'.format(CA_CO2))
                f.write('        thetaProperties   ((water H2) {} 0 0 0);\n'.format(CA_water))
                f.write('        value            uniform    0;\n')
                f.write('    }\n')
        f.write('}')

H2(file_path = 'E:\\generated_openfoampar\\H2_Control_CA150\\', Total_patch=185, Inlet_patch=141, Outlet_patch=158, up=None, down=None, CA_H2=150, CA_water=30, compressible = True)                    
#VE()

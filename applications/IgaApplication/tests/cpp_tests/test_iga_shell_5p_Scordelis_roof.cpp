//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

#include "test_creation_utility_Scordelis_roof.h"
#include "test_creation_utility.h"

#include "custom_elements/shell_5p_element.h"
#include "custom_utilities/director_utilities.h"


namespace Kratos
{
namespace Testing
{
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    ///@}
    ///@name Operations
    ///@{

    typename Shell5pElement::Pointer GetShell5pElementScordelis(
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
    {
        // Set the element properties
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 4.32e8);
        p_elem_prop->SetValue(POISSON_RATIO, 0.0);
        p_elem_prop->SetValue(THICKNESS, 0.25);

        auto p_quadrature_point = TestCreationUtilityScordelisRoof::GetQuadraturePointGeometry(
            rModelPart, PolynomialDegree, IntegrationPoint);

        return Kratos::make_intrusive<Shell5pElement>(1, p_quadrature_point, p_elem_prop);
    }

    Parameters GetDirectorParametersScordelisLoTest()
    {
        return Parameters(R"(
        {
            "model_part_name": "ModelPart",
            "brep_ids" : [1] ,
            "linear_solver_settings" : {
                "solver_type": "skyline_lu_factorization"
            }
        })");
    }

    // Tests the stiffness matrix of the Shell5pElement with a polynomial degree of p=4 (Scordelis Lo Roof).
    KRATOS_TEST_CASE_IN_SUITE(IgaShell5pElementP4Scordelis, KratosIgaFast5PSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DIRECTORINC);

        IntegrationPoint<3> integration_point(0.953089922969332, 0.953089922969332, 0.0, 0.014033587215607);

        auto p_shell_5p_element = GetShell5pElementScordelis(r_model_part, 4, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);
        TestCreationUtility::AddDirectorInc2DDofs(r_model_part);

        DirectorUtilities(r_model_part, GetDirectorParametersScordelisLoTest()).ComputeDirectors();

        p_shell_5p_element->Check(r_process_info);
        p_shell_5p_element->Initialize(r_process_info);

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_5p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        //const double tolerance = 1.0e-8;

        // const std::array<double, 125> expected_LHS_row_72{ 3.31347494115561e-05,-0.0003960152217280509,0.001084875366429851,1.715367327097879e-05,-5.953618914724417e-05,0.002692845355476361,-0.03218397217025212,0.08894284619206648,0.001394070883813759,-0.004838500257169575,0.08206734890058041,-0.9808411468399856,2.734264000841506,0.04248580460865035,-0.1474591656360989,1.11159574790744,-13.28541700812595,37.35552942585122,0.5754668620600415,-1.997331989839015,5.646178672087411,-67.48122733040603,191.3677348902012,2.922994918946871,-10.14518594099407,0.001757274768188859,-0.02100234585606551,0.0568541332798362,0.001230866283936651,-0.009934502952959339,0.1428128862252028,-1.706850899602881,4.676163353443921,0.1000319185982758,-0.8073760822662783,4.352375801123312,-52.01811764458926,144.2072346515098,3.048579952095885,-24.60576568176297,58.95258587778868,-704.5813795102861,1976.249715928131,41.29277425557886,-333.2846946063676,299.4405418284202,-3578.812605334494,10154.7483169131,209.7402601189065,-1692.875905991343,0.03311263023854576,-0.3957507027313587,1.044454782080836,0.03585174567335821,-0.4285596314502323,2.691047740785871,-32.16247820701539,86.50358950682951,2.913654351903946,-34.82899932175045,82.01256466587009,-980.1861979751957,2685.698805486688,88.79673977085103,-1061.456012804798,1110.853699782672,-13276.54716871586,37046.82133166722,1202.744814833775,-14377.40438569021,5642.409553409494,-67436.18141912093,191572.9177895399,6109.156259603934,-73028.1404875567,0.2130604069750364,-2.546421940575445,6.132092163462885,0.5079899221076598,-7.265643402876953,17.31531813421629,-206.9465938634658,521.3263510627027,41.28409982479887,-590.4781291408725,527.7028819183344,-6306.922081912083,16588.14069977122,1258.17719813656,-17995.53743830526,7147.693785131813,-85426.79840394152,234168.682736282,17041.91060511694,-243749.2614163977,36305.60505489032,-433912.4015703967,1237594.233870779,86561.74906333575,-1238092.415483469,-0.2479634467311842,2.963600316167391,-14.5901857236078,2.921766972672136,-43.36173247278853,-20.15187160658294,240.8498642464762,-1053.585018359659,237.4506148943586,-3524.003762992964,-614.14988973423,7340.141751421348,-28081.6453025252,7236.562032632952,-107398.2906106908,-8318.611666540191,99421.42190578765,-325811.6635167199,98018.66019442663,-1454708.098022836,-42253.10132880032,504994.6309862123,-1377821.503538731,497870.6240089242,-7388999.064980459 };
        // const std::array<double, 125> expected_LHS_row_73{ -9.807261086625984e-20,-0.000257765153951399,-0.0003484445150380983,0.000258083736117855,6.57358759085676e-07,-5.865645460722085e-18,-0.01545357362064315,-0.02088999573179635,0.02098242124037044,3.941006696589485e-05,-1.146197872543919e-16,-0.3035019222190682,-0.4102710489746547,0.6397079817025317,0.0007739977414006965,-6.837170109383887e-16,-1.84264592143276,-2.490871456588063,8.668128780557435,0.004699158973887663,9.400950271963094e-16,2.16185918242629,2.922380945809376,44.04535795121599,-0.005513224140996766,-7.050280345370785e-18,-0.01849794835788827,-0.02500535291895758,0.01844002990395988,0.0001097167121912405,-4.219359982071507e-16,-1.108991663137601,-1.499124518266478,1.499192213888682,0.006577752125421376,-8.255939239560082e-15,-21.78014676408387,-29.44219790895681,45.70715794243296,0.1291843856282504,-4.947860211610766e-14,-132.2334247823243,-178.7519021240216,619.3391512689452,0.7843149049950477,6.536528023404589e-14,155.1410611578966,209.7182299041545,3147.052980521583,-0.9201867594586216,-2.057355798547199e-16,-0.5388509003345507,-0.7284135879759632,0.534810548558616,0.00473405656084179,-1.232028063559415e-14,-32.30526675626565,-43.66995629745122,43.48068732536072,0.2838168405068799,-2.413865924636739e-13,-634.462345022194,-857.6602412398851,1325.635234415661,5.574047710034987,-1.453362227656772e-12,-3852.000158054865,-5207.097648478799,17962.59563527973,33.84161854288459,1.84297869735092e-12,4519.306620733526,6109.156259603934,91273.63573931654,-39.70421714984927,-2.920494991946382e-15,-7.635878791859763,-10.32210926004738,7.545278306073338,0.0802767788162583,-1.749997999288642e-13,-457.787304681848,-618.8325804524396,613.4404993279347,4.812765001195426,-3.433201378902376e-12,-8990.75711218208,-12153.62105294611,18702.55014195508,94.52075389868762,-2.076566751525869e-11,-54585.42668272418,-73788.06730484474,253423.3494911718,573.8622028773497,2.524996415729348e-11,64041.6069783791,86570.8430475022,1287726.482903803,-673.2759985532477,-1.682859679449362e-14,-43.92331737810367,-59.37512804994385,43.21034896206039,0.4791986550202635,-1.009015360076747e-12,-2633.297047177217,-3559.666225202514,3513.061131971741,28.72898675693875,-1.982098647229486e-11,-51716.88667043849,-69910.40184802846,107106.2568207621,564.2256553842019,-1.204298988699791e-10,-313987.8311007808,-424445.8021521573,1451314.023002052,3425.573370544011,1.402605129228078e-10,368381.9381357756,497975.2453534395,7374611.43603309,-4019.007211326162 };
        // const std::array<double, 125> expected_LHS_row_74{ -1.21030844949758e-19,-0.0003170726770943169,-0.0004288498657569198,-2.331608993271993e-07,4.938885646311937e-05,-9.835724726596162e-18,-0.02576824156731329,-0.03485250227467667,-1.894887560592176e-05,0.004014566431154071,-2.997422292902915e-16,-0.7853116037293411,-1.062169722837997,-0.0005774872970180123,0.1223713015230601,-4.059821896323422e-15,-10.63694045126246,-14.38703675298313,-0.007822019748861641,1.657822295136609,-2.062041004661029e-14,-54.02847495083736,-73.07688995465755,-0.03973073949027235,8.422240234028841,-6.419120126108898e-18,-0.01681578595755779,-0.02274363961317865,-1.236550064791084e-05,0.008160030370371342,-5.216505286823556e-16,-1.366603451328577,-1.848370339358155,-0.001004938367696893,0.6632882273421102,-1.589699607929516e-14,-41.64847686514933,-56.33123649339838,-0.03062657403557249,20.21827293723746,-2.153118651111916e-13,-564.1222317997975,-763.0044302698243,-0.4148345222169443,273.9070982469144,-1.093583944014884e-12,-2865.35591774109,-3875.574615744285,-2.107087793563751,1391.533836352301,-1.209702477473284e-16,-0.3168663575140581,-0.4285596314502274,-0.0002330052522686165,0.3488286703023415,-9.830369732765897e-15,-25.75136024462323,-34.82899932175265,-0.01893622624344457,28.3545957648475,-2.995659166261586e-13,-784.7936582830406,-1061.456012804701,-0.5771017941328249,864.3029583807636,-4.057256134097588e-12,-10629.87793328789,-14377.40438569081,-7.816798142824679,11709.16712400141,-2.060647545512223e-11,-53992.36331545447,-73028.1404875567,-39.70421714984927,59486.3417120566,-7.786721395928442e-16,-2.038928856835114,-2.757471936144805,-0.001499252506316627,5.860994330237037,-6.327062609116121e-14,-165.6996880854848,-224.1007240120327,-0.1218435394862892,476.4127344300074,-1.927883918878905e-12,-5049.781383376177,-6829.779549047639,-3.713312480450074,14522.00805321627,-2.610820995731597e-11,-68397.61285584453,-92509.76842244374,-50.29652375964095,196737.6368417132,-1.325881055164562e-10,-347408.6119462521,-469893.8591953657,-255.4734131262585,999492.3605779923,9.024400170485263e-16,2.371938026690263,3.209936851643524,0.001744856420134489,34.66474346106089,7.340810990377634e-14,192.7840646250735,260.856878780954,0.1418036529715935,2817.739453161042,2.23923999526816e-12,5875.843114170854,7949.491786634397,4.321618335975236,85890.45953007591,3.03580841950337e-11,79595.17257044488,107669.8026904518,58.5359784435616,1163607.427979147,1.543401736022676e-10,404328.6631067386,546864.5052893611,297.3244488100829,5911521.788459131 };
        // const std::array<double, 125> expected_RHS{-3.12052e-17,-1.63837e-17,1.21266e-17,0,0,-2.53603e-15,-9.82241e-16,7.27017e-16,0,0,-7.72883e-14,-1.92908e-14,1.42783e-14,0,0,-1.04686e-12,-1.1712e-13,8.66877e-14,0,0,-5.31738e-12,1.37409e-13,-1.01705e-13,0,0,-1.65494e-15,-1.17574e-15,8.7024e-16,0,0,-1.34496e-13,-7.04883e-14,5.21728e-14,0,0,-4.09892e-12,-1.38436e-12,1.02465e-12,0,0,-5.55196e-11,-8.40485e-12,6.22095e-12,0,0,-2.82003e-10,9.86088e-12,-7.29865e-12,0,0,-3.11844e-14,-3.42498e-14,2.53504e-14,0,0,-2.53434e-12,-2.05335e-12,1.51981e-12,0,0,-7.72367e-11,-4.03269e-11,2.98484e-11,0,0,-1.04616e-09,-2.44836e-10,1.81218e-10,0,0,-5.31383e-09,2.8725e-10,-2.12612e-10,0,0,-2.00653e-13,-4.85342e-13,3.59232e-13,0,0,-1.6307e-11,-2.90973e-11,2.15367e-11,0,0,-4.96973e-10,-5.71459e-10,4.22972e-10,0,0,-6.73146e-09,-3.46949e-09,2.56798e-09,0,0,-3.41914e-08,4.07053e-09,-3.01285e-09,0,0,2.33524e-13,-2.7918e-12,2.06638e-12,0,0,1.89784e-11,-1.67374e-10,1.23884e-10,0,0,5.78386e-10,-3.28716e-09,2.43303e-09,0,0,7.83419e-09,-1.99573e-08,1.47716e-08,0,0,3.97926e-08,2.34146e-08,-1.73306e-08,0,0 };

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          // KRATOS_CHECK_NEAR(left_hand_side_matrix(72,i), expected_LHS_row_72[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          // KRATOS_CHECK_NEAR(left_hand_side_matrix(73,i), expected_LHS_row_73[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          // KRATOS_CHECK_NEAR(left_hand_side_matrix(74,i), expected_LHS_row_74[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          // KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell5pElement with a polynomial degree of p=5 (Scordelis Lo Roof).
    KRATOS_TEST_CASE_IN_SUITE(IgaShell5pElementP5Scordelis, KratosIgaFast5PSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DIRECTORINC);

        IntegrationPoint<3> integration_point(0.619309593041599, 0.966234757101576, 0.0, 0.020041279329452);
        auto p_shell_5p_element = GetShell5pElementScordelis(r_model_part, 5, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);
        TestCreationUtility::AddDirectorInc2DDofs(r_model_part);

        DirectorUtilities(r_model_part, GetDirectorParametersScordelisLoTest()).ComputeDirectors();

        p_shell_5p_element->Initialize(r_process_info);

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_5p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        //const double tolerance = 1.0e-8;

        // const std::array<double, 180> expected_LHS_row_105{ 0.01535086990253638,0.0007026666752406215,-0.0005390487579320884,-2.55982033322816e-21,4.665206191192949e-20,0.1248424827150979,0.003869743921223019,-0.002968663133962052,-2.081652444351261e-20,3.795872873371056e-19,0.4061175713452437,0.006585229682705845,-0.005051840376440981,-6.771219962634575e-20,1.235413732173715e-18,0.6605583013564119,0.0009432401813355855,-0.0007236040445585216,-1.10127454547301e-19,2.010403131756449e-18,0.537205586185999,-0.00717943299528232,0.005507681771642725,-8.955591143809286e-20,1.635776038932892e-18,0.1747551017239875,-0.004921447465223026,0.003775474541251055,-2.913083142614286e-20,5.323833711879475e-19,1.5802397579716,0.09112982861863964,-0.06990999097917933,-2.358657734626097e-19,1.310146003013376e-17,12.85086957032281,0.5018725275655275,-0.3850100939137235,-1.917873642793374e-18,1.066007078230582e-16,41.80247706103692,0.8540476922345305,-0.6551802781256709,-6.237850981095293e-18,3.469448379866791e-16,67.98944588210037,0.1223301447188924,-0.09384522547023057,-1.014425019910648e-17,5.645868189903064e-16,55.29055750510337,-0.9311107549180946,0.7142989892952166,-8.248497758462385e-18,4.59378856172215e-16,17.98541264087753,-0.6382694382195239,0.4896465991936018,-2.682808712788876e-18,1.495103365468038e-16,63.48380412804771,4.946349402836643,-3.794577992451629,-7.58507799117205e-18,1.063311779579898e-15,516.2245365405411,27.2406621920949,-20.89759716421438,-6.166111197690955e-17,8.65168351919225e-15,1679.09123806361,46.35604342192551,-35.56190795682458,-2.005037484099958e-16,2.815792019454611e-14,2730.737239936885,6.639841723080576,-5.093735849172941,-3.259894121183784e-16,4.582162491252961e-14,2220.524124460706,-50.53887619867481,38.77075632318847,-2.650051708137274e-16,3.728295795864264e-14,722.2557881076956,-34.64402054126398,26.57706263947561,-8.617177578913038e-17,1.213417310743478e-14,1178.854455390904,141.5459894118657,-108.5866065252154,-6.775818691417998e-17,3.823489903608274e-14,9584.411255712042,779.5257003080372,-598.0109422199287,-5.501126313017697e-16,3.110996788654165e-13,31169.55416696892,1326.536299197334,-1017.648580230853,-1.786487215423643e-15,1.012509549710526e-12,50683.38959089881,190.0073952908913,-145.7636373525096,-2.900792544467523e-15,1.647664013774221e-12,41206.96461024128,-1446.233303131143,1109.473799174977,-2.355059505408358e-15,1.340627575329893e-12,13400.94719434576,-991.3820810770063,760.5359671535399,-7.647958752025828e-16,4.363224809533756e-13,8162.644928045262,2135.496608019622,-1638.240199347676,1.230191671725339e-15,6.458848006910019e-13,66328.40055930353,11760.66164635808,-9022.158409259897,1.002202789631464e-14,5.255261120435929e-12,215589.5445813974,20013.38066250671,-15353.20852444297,3.265858542517027e-14,1.71038344987692e-11,350368.8019215446,2866.631190529059,-2199.128032042965,5.321181669399049e-14,2.783316759135485e-11,284703.4454515417,-21819.24282047274,16738.57058521637,4.33498660895464e-14,2.264653418193691e-11,92537.88426245023,-14956.92728694087,11474.16457987722,1.412623353504472e-14,7.370565502699912e-12,-9218.792199581252,13483.82925610458,-10344.08158061931,2.198993096936098e-14,4.190900750603706e-12,-75407.82900545046,74258.49003121104,-56967.19117002255,1.789653630935572e-13,3.409936051212728e-11,-246720.5061668415,126367.3314568834,-96942.34188877542,5.826046751618624e-13,1.109800825690517e-10,-403599.4992817534,18100.31697927344,-13885.60711593906,9.483066969490291e-13,1.805983658981607e-10,-330105.4542768003,-137769.8019204889,105689.7149425197,7.717800402815789e-13,1.469442352032354e-10,-107994.4966109816,-94440.16580298342,72449.50681283661,2.512454050712927e-13,4.782458829446752e-11 };
        // const std::array<double, 180> expected_LHS_row_106{ 5.75202005418484e-05,0.02500523939263251,-0.007409532607479236,-9.684625105296157e-06,-0.000443241940292455,0.0004678711538052481,0.203382425767615,-0.06026937554141768,-7.87750508431073e-05,-0.003605344534604896,0.001522271581120741,0.661691901297655,-0.1960932123835109,-0.0002563035148123079,-0.01173039645521905,0.002476441161034286,1.07638644123565,-0.3190056939984365,-0.0004169561999792394,-0.01908308619037302,0.002014345173398602,0.8754888530587325,-0.2594802532083494,-0.0003391535087368576,-0.01552224513176344,0.0006553899267123456,0.2848347768893942,-0.08442482760874999,-0.0001103474201794664,-0.005050338116797302,0.005920369568736848,2.573821012050619,-0.762639402863264,-0.001256333370820218,-0.1241567917852117,0.04815647572478675,20.93411693298212,-6.203333327486725,-0.01021905588353203,-1.009895431756284,0.1566825264745599,68.10687265035553,-20.18324478249228,-0.03324885116541202,-3.28580907558459,0.2548921379001088,110.7891509239748,-32.83423189150507,-0.05408944409371538,-5.345375838549674,0.2073300814875393,90.11011244436683,-26.70746937678482,-0.04399652709584889,-4.34794630473167,0.06745718097664608,29.31634463968616,-8.689576450197899,-0.0143147664311178,-1.414653535272884,0.2377847280764802,103.3821745519335,-30.63052079896628,-0.06820873425892862,-10.05243323403156,1.934148595352496,840.8365392001872,-249.1496370832197,-0.554812029455459,-81.7668228588477,6.292970653188146,2735.50477238712,-810.6364506971761,-1.805143528152044,-266.0376115775432,10.23744497631784,4449.72711550134,-1318.748569919496,-2.936619056789881,-432.7917381481994,8.327170538292508,3619.089141160335,-1072.674309901208,-2.388655347917109,-352.0342245853362,2.709338876417475,1177.403086589127,-349.0066875120175,-0.7771759647220539,-114.5383188795411,4.413294227835287,1919.073601449422,-568.5037134835208,-1.952376613821966,-360.6205913103415,35.89787662435336,15607.59181414816,-4624.227407969473,-15.88069391913896,-2933.299761419232,116.7977913649668,50773.85612512227,-15045.44551843875,-51.66962922371957,-9543.822729704101,190.0073952908913,82587.51396777207,-24476.02711516594,-84.05648385708594,-15525.95365357471,154.5526239978969,67167.41495070975,-19908.87882343361,-68.37181323555691,-12628.86181358611,50.28543978104911,21850.57437149278,-6477.57831189142,-22.2455407629711,-4108.943109456506,30.5070157591473,13272.53903460654,-3929.7973290284,-29.46291514954061,-6077.64128904402,248.1450434443923,107925.8501202788,-31965.09722670799,-239.6522956389005,-49435.73431948997,807.3678929746681,351040.2336056178,-104002.049486406,-779.7357798935947,-160844.7561467038,1313.431260923873,570897.6749216857,-169191.1997116551,-1268.47916232901,-261663.3082929716,1068.349195071473,464226.1666360231,-137620.663252462,-1031.785014075864,-212837.7963546167,347.5994630452207,150994.4233924553,-44776.43514869614,-335.7028942645962,-69249.18568456934,-35.16407260482833,-15203.69985499062,4529.701234377816,-186.0804921450748,-39344.37353438357,-286.0256930109767,-123878.3217312494,36844.7357934564,-1513.584683997098,-320028.4295174995,-930.6168597908789,-403738.4122443792,119878.5072522402,-4924.618522275967,-1041248.712321579,-1513.933469770143,-657920.7411437159,195019.1281270858,-8011.400963012826,-1693910.260638206,-1231.438334034323,-536063.0075256441,158629.1871962643,-6516.499206975617,-1377832.182235363,-400.6623542735907,-174709.6300790138,51611.79679595254,-2120.216531944235,-448293.2931023353 };
        // const std::array<double, 180> expected_LHS_row_107{ -4.412645960116244e-05,-0.00740953260747924,0.02103088723746291,-1.264902112031777e-05,-0.0005801443369644487,-0.000358926035939601,-0.0602693755414177,0.1710549180470909,-0.0001028875429905555,-0.00471891231933648,-0.001167806349657129,-0.196093212383511,0.556510708024565,-0.0003347562282302113,-0.0153535126060505,-0.001899794851506678,-0.3190056939984366,0.9052770022712685,-0.00054458357679748,-0.02497719107679893,-0.001545299218004319,-0.2594802532083496,0.7363078873765531,-0.0004429660258812048,-0.02031652593386472,-0.0005027805336499228,-0.08442482760875007,0.2395506795877686,-0.0001441239937784837,-0.006610210487755366,-0.00454179481534928,-0.7626394028632643,2.164753704441047,-0.001640223278202255,-0.1623581897942553,-0.03694310451957051,-6.203333327486725,17.60675035193136,-0.01334162868767539,-1.320626632243325,-0.1201985582379466,-20.18324478249228,57.28091011681543,-0.04340849405232854,-4.296807457745278,-0.1955397846278378,-32.83423189150509,93.17740573108381,-0.07061721623267032,-6.990073453540418,-0.1590526872069842,-26.70746937678483,75.78466250038424,-0.0574402698987088,-5.685747775334073,-0.05174958611294547,-8.689576450197903,24.65539789420659,-0.01868883981567623,-1.849922063273083,-0.1824158834357051,-30.63052079896629,86.95246396148988,-0.08901476628435273,-13.13438783933306,-1.483776639362954,-249.1496370832198,707.1967471164201,-0.7240489604490841,-106.8355245934255,-4.827634686307789,-810.6364506971765,2300.692633003957,-2.355774975360014,-347.6014095082169,-7.853627037304065,-1318.748569919497,3742.371942288943,-3.832389822892967,-565.4801638093502,-6.388165390393934,-1072.674309901209,3043.724178052184,-3.117278158564815,-459.9633472515902,-2.078461676939291,-349.0066875120175,990.2016051782931,-1.0142416160214,-149.6542544605763,-3.385646218512977,-568.5037134835208,1614.137503776478,-2.546886965870179,-470.793031117514,-27.53895479692531,-4624.227407969474,13127.23183215597,-20.71646016720682,-3829.445351266828,-89.60109619961327,-15045.44551843876,42703.72512099721,-67.40334025187035,-12459.53166602582,-145.7636373525096,-24476.02711516595,69458.97317196564,-109.652185760104,-20269.244642016,-118.5646095606899,-19908.87882343361,56488.61783210951,-89.19143915664091,-16487.0674665589,-38.57633329026368,-6477.578311891421,18376.10726432881,-29.01944093001275,-5364.253053281149,-23.40337108992944,-3929.797329028401,11164.65999399898,-38.41897872085,-7927.890173636464,-190.363770147346,-31965.097226708,90780.29455245027,-312.5012036257156,-64485.70850370979,-619.3700017909825,-104002.049486406,295255.2345489805,-1016.757920374505,-209811.5140422208,-1007.595087084143,-169191.1997116551,480146.2814355438,-1654.068299012827,-341322.693400367,-819.5810715568975,-137620.663252462,390408.692772531,-1345.424453048076,-277632.954421186,-266.6599476177311,-44776.43514869614,126977.0756107071,-437.7490240126621,-90331.00801199953,26.97601911315307,4529.701234377818,-12774.04216714242,-242.5464890384625,-51280.1393735537,219.4238036141898,36844.7357934564,-104115.4078114744,-1972.880911555584,-417114.2696458157,713.9200990414914,119878.5072522402,-339437.5408000004,-6418.990613484531,-1357128.245496903,1161.409791053436,195019.1281270859,-553315.6691059987,-10442.45505511228,-2207784.783159387,944.6944444944063,158629.1871962644,-450976.9048786437,-8493.92639310953,-1795819.099972325,307.3669949515807,51611.79679595254,-147025.9060582684,-2763.594775015166,-584289.9660192849 };
        // const std::array<double, 180> expected_RHS{ -1.918281691746684e-16,-8.785481996251462e-18,6.739757732629536e-18,0,0,-1.560336473076147e-15,-4.838363159654087e-17,3.711736650585826e-17,0,0,-5.076730742281859e-15,-8.233550680168764e-17,6.316345180303315e-17,0,0,-8.258858097067607e-15,-1.179338642810721e-17,9.047263133273145e-18,0,0,-6.717781632519783e-15,8.976486511438141e-17,-6.886286307697602e-17,0,0,-2.18570603983345e-15,6.153314170820405e-17,-4.720497609781776e-17,0,0,-1.974425757403166e-14,-1.13940150694691e-15,8.740886521982217e-16,0,0,-1.606004235932615e-13,-6.274941178660958e-15,4.813803425779075e-15,0,0,-5.225315960678757e-13,-1.067820758896416e-14,8.191756832357246e-15,0,0,-8.500577482309266e-13,-1.529500859932995e-15,1.173352270506613e-15,0,0,-6.914396954797441e-13,1.164173150953492e-14,-8.930921490256671e-15,0,0,-2.249677053024137e-13,7.980319624969928e-15,-6.122079690584455e-15,0,0,-7.930050419664578e-13,-6.184449207144953e-14,4.744382765077011e-14,0,0,-6.450328414420075e-12,-3.405915716333977e-13,2.612838634879761e-13,0,0,-2.098687118089465e-11,-5.795922864298589e-13,4.446325877060351e-13,0,0,-3.41415764954555e-11,-8.301832429453497e-14,6.368727331667865e-14,0,0,-2.777086769027732e-11,6.318904860575413e-13,-4.847530040353483e-13,0,0,-9.035565096105557e-12,4.331561883716927e-13,-3.322945481261105e-13,0,0,-1.471820584385552e-11,-1.769757675207382e-12,1.357664608663325e-12,0,0,-1.197183578158391e-10,-9.746454822731354e-12,7.476965325892853e-12,0,0,-3.895171830092248e-10,-1.658575991238141e-11,1.272371893411254e-11,0,0,-6.336690488723375e-10,-2.375673429953414e-12,1.822491170838721e-12,0,0,-5.154284342435074e-10,1.80823384611015e-11,-1.387181494643465e-11,0,0,-1.677004558133223e-10,1.239530737917193e-11,-9.509025093072761e-12,0,0,-1.017399961219246e-10,-2.670023734423969e-11,2.048301176653311e-11,0,0,-8.275563875872514e-10,-1.470442313529407e-10,1.128045673216793e-10,0,0,-2.692548066605946e-09,-2.502284535329203e-10,1.91962052320181e-10,0,0,-4.380254445331708e-09,-3.58416552281524e-11,2.749582471141023e-11,0,0,-3.562912997507528e-09,2.728072523226714e-10,-2.092833141247857e-10,0,0,-1.159233938232715e-09,1.870073251355808e-10,-1.434621419950176e-10,0,0,1.172711430934898e-10,-1.685890954334341e-10,1.293326490304477e-10,0,0,9.538872345775769e-10,-9.284581868275071e-10,7.122640791664981e-10,0,0,3.103579729122876e-09,-1.579978038731676e-09,1.212075695843459e-09,0,0,5.048923387305828e-09,-2.263093078854139e-10,1.736125472043029e-10,0,0,4.106810456918425e-09,1.72254378505299e-09,-1.321444606005378e-09,0,0,1.336197112553485e-09,1.180790843825042e-09,-9.05840365239331e-10,0,0 };

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          // KRATOS_CHECK_NEAR(left_hand_side_matrix(105,i), expected_LHS_row_105[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          // KRATOS_CHECK_NEAR(left_hand_side_matrix(106,i), expected_LHS_row_106[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          // KRATOS_CHECK_NEAR(left_hand_side_matrix(107,i), expected_LHS_row_107[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          // KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }
}
}

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

#include "custom_elements/shell_3p_element.h"

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

    typename Element::Pointer GetShell3pElementScordelis(
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
    {
        // Set the element properties
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 4.32e8);
        p_elem_prop->SetValue(POISSON_RATIO, 0.0);
        p_elem_prop->SetValue(THICKNESS, 0.25);
        const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStress2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        auto p_quadrature_point = TestCreationUtilityScordelisRoof::GetQuadraturePointGeometry(
            rModelPart, PolynomialDegree, IntegrationPoint);

        return Kratos::make_intrusive<Shell3pElement>(1, p_quadrature_point, p_elem_prop);
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=4 (Scordelis Lo Roof).
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP4Scordelis, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.953089922969332, 0.953089922969332, 0.0, 0.014033587215607);
        auto p_shell_3p_element = GetShell3pElementScordelis(r_model_part, 4, integration_point);

        TestCreationUtilityScordelisRoof::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize(r_process_info);

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

        //Check RHS and LHS results
        const double tolerance = 1.0e-6;

        const std::array<double, 75> expected_LHS_row_72{-0.0121324873940569,-0.00325599890435487,0.00240996641364744,-0.859568758238585,-0.195204115086357,0.144482653398527,-22.3431062826556,-3.833729699559,2.83758075059011,-250.444682486708,-23.2756561906353,17.2277544688395,-1006.99930056829,27.3078460041853,-20.2122278392417,-0.754515464569957,-0.233659587662703,0.172945929969653,-52.2459471150918,-14.0083932399781,10.3684792929014,-1315.73676928285,-275.119165307797,203.632730728415,-14076.183742387,-1670.32618494086,1236.31184278328,-52473.7873098818,1959.68740307631,-1450.48599873456,-18.5959049070138,-6.80657534272742,5.03796790164583,-1246.97624504774,-408.068784903039,302.036976965475,-29947.9864461454,-8014.30553579081,5931.88380503813,-296539.4280861,-48657.113275156,36014.1425788958,-952050.897379708,57086.2941711929,-42253.1013288009,-215.571944827077,-96.4537393784795,71.3913853170952,-13774.0629232677,-5782.60846985042,4280.06661088661,-305635.268027351,-113568.086523072,84058.7734279465,-2593733.19254302,-689503.940886444,510344.563504426,-5321431.3012989,808951.089618747,-598754.794928575,-964.231291487823,-554.823920403535,410.659540417027,-56818.4153545866,-33262.883555108,24619.9198905316,-1075019.40509719,-653269.550807246,483525.247627791,-5667691.91082686,-3966184.01817668,2935618.39390265,16384286.6104441,4653271.27645944,-3444174.22096139};
        const std::array<double, 75> expected_LHS_row_73{-0.00364908072639826,-0.00991775088125891,0.0073642940391869,-0.296559058400869,-0.765346303345953,0.568023991409893,-9.03795521191612,-22.0850516162108,16.3831145833325,-122.418388286365,-282.343254631751,209.345280317556,-621.805268970262,-1348.77349219962,999.563813108502,-0.193526059543626,-0.561917431635313,0.417042645690037,-15.7277709915,-42.7459439143073,31.7098703147866,-479.32068091579,-1213.68466420156,899.90219307765,-6492.36070042372,-15232.7826676977,11289.0262771341,-32976.9420108317,-71242.4593569498,52771.7258679007,-3.64664477473573,-12.0034089951456,8.90428258830844,-296.361089765612,-890.350870038262,660.158577977399,-9031.92190553705,-24538.2039523195,18185.2980387821,-122336.6676288,-297191.471158319,220141.341021844,-621390.181529813,-1330806.67828469,985289.18905296,-23.4640260891792,-108.214744698997,80.2337737516887,-1906.91026179862,-7586.70098479888,5622.34724504934,-58115.1344093582,-194390.167604492,143989.052724168,-787164.897659499,-2134078.66225517,1579975.16518729,-3998282.34764983,-8304716.21870224,6145230.08517709,27.3078460041855,-266.399586960411,197.392744571164,2219.29568161417,-14696.0764388357,10883.6697864567,67635.4149510236,-235877.119357513,174577.515968701,916116.344377013,-322469.150181996,237789.612251406,4653271.27645944,12957013.6291438,-9588848.61468};
        const std::array<double, 75> expected_LHS_row_74{0.00270091061134752,0.0073642940391869,-0.00541894935257825,0.21950172325092,0.568023991409894,-0.418343971922984,6.68955031883932,16.3831145833325,-12.0767085842989,90.6094298092045,209.345280317556,-154.4555219987,460.236584245489,999.563813108502,-738.14616166423,0.14324062057941,0.417042645690037,-0.307148667382859,11.6410972375813,31.7098703147866,-23.3745809177562,354.774917408195,899.90219307765,-663.939998882923,4805.39818744374,11289.0262771341,-8336.38657017108,24408.2768469635,52771.7258679007,-39004.5395044392,2.69910761267817,8.90428258830844,-6.56383801835856,219.355194404953,660.158577977399,-487.064091431539,6685.08469518161,18185.2980387821,-13428.9182311019,90548.943290069,220141.341021844,-162708.492037047,459929.35232693,985289.189052959,-728899.627002229,17.3671786953728,80.2337737516887,-59.2004302127478,1411.42236830232,5622.34724504934,-4152.04394283289,43014.6096989887,143989.052724168,-106428.148833675,582629.485170998,1579975.16518729,-1168881.63853518,2959376.35520325,6145230.08517709,-4550633.29450364,-20.2122278392418,197.392744571165,-145.813584242182,-1642.63816166811,10883.6697864567,-8047.31010598395,-50061.1588618974,174577.515968701,-129228.798210361,-678074.436078321,237789.612251408,-177204.951684685,-3444174.22096139,-9588848.61468,7099245.51498888};
        const std::array<double, 75> expected_RHS{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(72,i), expected_LHS_row_72[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(73,i), expected_LHS_row_73[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(74,i), expected_LHS_row_74[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=5 (Scordelis Lo Roof).
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP5Scordelis, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.619309593041599, 0.966234757101576, 0.0, 0.020041279329452);
        auto p_shell_3p_element = GetShell3pElementScordelis(r_model_part, 5, integration_point);

        TestCreationUtilityScordelisRoof::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize(r_process_info);

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

        //Check RHS and LHS results
        const double tolerance = 1.0e-6;

        const std::array<double, 108> expected_LHS_row_105{-0.0365121097495112,-0.00148169014127313,0.00113667441254417,-0.286063463442857,-0.00816000191750132,0.00625992245447856,-0.895186849811542,-0.0138860575614221,0.0106526498904973,-1.39845801648961,-0.00198897959271944,0.00152584008433396,-1.09046487130985,0.0151390345719724,-0.0116138676699276,-0.339488309481346,0.0103776946409439,-0.0079612191719262,-3.86941189807431,-0.192162476744624,0.147416902011002,-30.0567555135908,-1.05828211650285,0.811858140628705,-93.1821392632363,-1.80090232018991,1.38155713521345,-144.088250363733,-0.257953558632089,0.197888344907188,-111.10046567761,1.96340266958451,-1.50621881989442,-34.1626374675073,1.345897802485,-1.0325017028659,-163.025044540636,-10.4302045389666,8.00150199208212,-1249.12999165348,-57.4414897333175,44.0660768264134,-3813.91438097237,-97.7494663499944,74.9882273930885,-5797.35793948495,-14.0012161774049,10.7409934979267,-4384.4185783252,106.569668454756,-81.7546205590273,-1318.76984498132,73.0527083449288,-56.0421791504824,-3320.22819989119,-298.4733792541,228.973011004997,-24805.6533187437,-1643.76024324751,1261.00603413033,-73546.1740384414,-2797.22352833731,2145.88213976748,-107994.496610979,-400.662354273576,307.36699495158,-78366.1978295033,3049.62466947285,-2339.51096325359,-22413.9194334742,2090.49483563967,-1603.71721660077,-29801.0978102749,-4503.05156387464,3454.50330568124,-209193.862761672,-24799.3209729072,19024.7291340905,-572585.883887545,-42201.5585284838,32374.8065891804,-755710.035847442,-6044.77104626932,4637.22905794497,-471719.919047284,46009.5207533043,-35296.0740690118,-106960.202754276,31539.1813582309,-24195.1940178852,-60048.0870747923,-28432.9079197312,21812.2251112593,-278745.317024063,-156586.43914968,120124.845117692,-224686.595295578,-266466.641718193,204419.132584221,744355.267500992,-38167.5439680175,29280.1236976731,1508238.35254338,290510.656707656,-222864.430863569,784451.17250335,199142.876047965,-152771.895647277};
        const std::array<double, 108> expected_LHS_row_106{-0.028589504465742,-0.0419772649372529,0.0323368591779443,-0.232547945158418,-0.338033942394766,0.260352732598042,-0.756620974136903,-1.08873369675989,0.838381251198945,-1.23087584823398,-1.75310242184669,1.34972463258284,-1.00119835793212,-1.41128989588255,1.0863556659016,-0.325751180629361,-0.454399134914215,0.349712486183971,-2.942625905852,-4.35707147816569,3.35333402822234,-23.9354133820609,-34.9975050163207,26.9300018879772,-77.8765676779785,-112.426609728337,86.4938982214083,-126.690099236978,-180.55061076956,138.877618324644,-103.050132557485,-144.952172431056,111.474324330799,-33.5285231729165,-46.5408625425934,35.7850184523907,-118.187132193356,-177.451535128073,136.445729249328,-961.337919258613,-1419.30794923465,1091.11965330959,-3127.82137227096,-4539.49322259554,3489.14059016127,-5088.35984769012,-7257.34172009897,5577.04064766197,-4138.88820012546,-5799.39367126299,4455.78002239364,-1346.63396818637,-1853.13693301312,1423.51923760015,-2193.55798218286,-3387.28275688987,2602.12937909089,-17842.4708953496,-26861.7866172858,20631.2064734532,-58052.491930867,-85150.9962279726,65386.9800658484,-94440.1658029814,-134868.492294533,103543.620322924,-76817.9333930741,-106728.509762083,81922.6558729156,-24993.5812399706,-33757.5848454581,25906.280900945,-15163.0289023081,-25578.5180948889,19631.2127254794,-123336.562822721,-197621.804769872,151638.219270497,-401289.421181739,-609029.85577898,467210.197559661,-652820.201353583,-935528.05814377,717509.95793418,-531006.043020412,-716018.514284338,549021.997900696,-172768.806565454,-218340.020314379,167374.621104699,17477.7452320946,-225.953396371241,170.828038060793,142164.539598656,64158.3232484073,-49249.5730962431,462548.367673529,423501.391161354,-324961.438008624,752476.64797934,1038366.60085574,-796559.616368056,612066.915944527,1128861.59427326,-865793.447486263,199142.876047965,459784.505147712,-352565.709528512};
        const std::array<double, 108> expected_LHS_row_107{0.0219323577098269,0.0323368591779443,-0.0246323023365684,0.178398500191258,0.260352732598041,-0.198384993031771,0.580439646144023,0.838381251198945,-0.639039699976212,0.94426293509924,1.34972463258284,-1.02913219503692,0.768066496254628,1.0863556659016,-0.828586477086793,0.24989909938882,0.349712486183971,-0.266819074634565,2.25742716354823,3.35333402822234,-2.55839792774278,18.3619848625562,26.9300018879772,-20.5526922714294,59.7427892313992,86.4938982214083,-66.0326967764886,97.1899522808589,138.877618324644,-106.058926754805,79.0546185228253,111.474324330799,-85.1591673091166,25.7213119797787,35.7850184523907,-27.3463636796714,90.6669250972625,136.445729249328,-104.264276655336,737.48767315871,1091.11965330959,-834.04920547777,2399.49934323942,3489.14059016126,-2667.9752525957,3903.5209046588,5577.04064766197,-4265.90829349792,3175.13640835948,4455.78002239364,-3409.38611301479,1033.06645030724,1423.51923760015,-1089.58448651052,1682.78181876602,2602.12937909088,-1991.54309911699,13687.8012199512,20631.2064734532,-15795.5447562606,44534.8054388382,65386.9800658484,-50078.4891062114,72449.506812835,103543.620322924,-79329.3871465127,58930.6609257745,81922.6558729156,-62786.5364725546,19173.7553500274,25906.280900945,-19861.8793016063,11632.2748527651,19631.2127254794,-15048.6565317794,94617.2962798237,151638.219270497,-116285.543025717,307848.047561424,467210.197559661,-358425.944040453,500809.176089244,717509.95793418,-550667.446990929,407359.788119404,549021.997900696,-421532.090351251,132539.102635863,167374.621104699,-128563.009894491,-13408.0029561496,170.828038060797,-134.324029554631,-109061.125556296,-49249.5730962431,37741.6574011364,-354842.675572379,-324961.438008624,249197.388477982,-577260.338021954,-796559.616368056,611105.045270188,-469545.408138557,-865793.447486263,664464.14355209,-152771.895647277,-352565.709528512,270674.022512257};
        const std::array<double, 108> expected_RHS{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(105,i), expected_LHS_row_105[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(106,i), expected_LHS_row_106[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(107,i), expected_LHS_row_107[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }
}
}

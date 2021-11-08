//
// Authors: Miquel Santasusana msantasusana@cimne.upc.edu
//          Miguel Angel Celigueta maceli@cimne.upc.edu
//

#include "spheric_continuum_particle.h"
#include <cmath>

namespace Kratos {

    SphericContinuumParticle::SphericContinuumParticle():SphericParticle() {
        mContinuumInitialNeighborsSize = 0;
        mInitialNeighborsSize = 0;
    }

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericParticle(NewId, pGeometry){
        mContinuumInitialNeighborsSize = 0;
        mInitialNeighborsSize = 0;
    }

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties){
        mContinuumInitialNeighborsSize = 0;
        mInitialNeighborsSize = 0;
    }

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes){
        mContinuumInitialNeighborsSize = 0;
        mInitialNeighborsSize = 0;
    }

    Element::Pointer SphericContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
      return SphericParticle::Pointer(new SphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Destructor

    SphericContinuumParticle::~SphericContinuumParticle() {
    }

    void SphericContinuumParticle::SetInitialSphereContacts(const ProcessInfo& r_process_info) {

        std::vector<SphericContinuumParticle*> ContinuumInitialNeighborsElements;
        std::vector<SphericContinuumParticle*> DiscontinuumInitialNeighborsElements;
        std::vector<int> DiscontinuumInitialNeighborsIds;
        std::vector<double> DiscontinuumInitialNeighborsDeltas;
        mIniNeighbourFailureId.clear(); // We will have to build this vector, we still don't know its size, it applies only to continuum particles
        size_t continuum_ini_size    = 0;
        size_t discontinuum_ini_size = 0;
        unsigned int neighbours_size = mNeighbourElements.size();
        mIniNeighbourIds.resize(neighbours_size);
        mIniNeighbourDelta.resize(neighbours_size);

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {

            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();

            double distance = DEM_MODULUS_3(other_to_me_vect);
            double radius_sum = GetRadius() + neighbour_iterator->GetRadius();
            double initial_delta = radius_sum - distance;
            int r_other_continuum_group = neighbour_iterator->mContinuumGroup; // finding out neighbor's Continuum Group Id
            if ((r_other_continuum_group  == this->mContinuumGroup) && (this->mContinuumGroup != 0)) {

                mIniNeighbourIds[continuum_ini_size]   = neighbour_iterator->Id();
                mIniNeighbourDelta[continuum_ini_size] = initial_delta;
                mIniNeighbourFailureId.push_back(0);
                array_1d<double, 3> vector_of_zeros(3,0.0);
                //mArrayOfOldDeltaDisplacements.push_back(vector_of_zeros);
                //mArrayOfDeltaDisplacements.push_back(vector_of_zeros);
                ContinuumInitialNeighborsElements.push_back(neighbour_iterator);
                continuum_ini_size++;
            } else {
                DiscontinuumInitialNeighborsIds.push_back(neighbour_iterator->Id());
                DiscontinuumInitialNeighborsDeltas.push_back(initial_delta);
                DiscontinuumInitialNeighborsElements.push_back(neighbour_iterator);
                discontinuum_ini_size++;
            }
        }

        mContinuumInitialNeighborsSize = continuum_ini_size;
        mInitialNeighborsSize = neighbours_size;

        for (unsigned int j = 0; j < continuum_ini_size; j++) {
            mNeighbourElements[j] = ContinuumInitialNeighborsElements[j];
        }

        for (unsigned int k = 0; k < discontinuum_ini_size; k++) {

            mIniNeighbourIds[continuum_ini_size + k]   = DiscontinuumInitialNeighborsIds[k];
            mIniNeighbourDelta[continuum_ini_size + k] = DiscontinuumInitialNeighborsDeltas[k];
            mNeighbourElements[continuum_ini_size + k] = DiscontinuumInitialNeighborsElements[k];
        }
    }//SetInitialSphereContacts

    void SphericContinuumParticle::CreateContinuumConstitutiveLaws() {

        mContinuumConstitutiveLawArray.resize(mContinuumInitialNeighborsSize);

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            Properties::Pointer properties_of_this_contact = GetProperties().pGetSubProperties(mNeighbourElements[i]->GetProperties().Id());
            mContinuumConstitutiveLawArray[i] = (*properties_of_this_contact)[DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER]-> Clone();
            SphericContinuumParticle* p_cont_neighbour_particle = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            mContinuumConstitutiveLawArray[i]->Initialize(this, p_cont_neighbour_particle, properties_of_this_contact);
        }
    }

    void SphericContinuumParticle::SetInitialFemContacts() {

        std::vector<DEMWall*>& rFemNeighbours = this->mNeighbourRigidFaces;

        unsigned int fem_neighbours_size = rFemNeighbours.size();

        mFemIniNeighbourIds.resize(fem_neighbours_size);
        mFemIniNeighbourDelta.resize(fem_neighbours_size);
        mContactConditionWeights.resize(fem_neighbours_size);

        for (unsigned int i = 0; i < rFemNeighbours.size(); i++) {

            double LocalCoordSystem[3][3]            = {{0.0}, {0.0}, {0.0}};
            array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
            array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);
            double DistPToB = 0.0;
            int ContactType = -1;
            array_1d<double, 4>& Weight = this->mContactConditionWeights[i];

            rFemNeighbours[i]->ComputeConditionRelativeData(i, this, LocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

            double initial_delta = -(DistPToB - GetRadius());

            mFemIniNeighbourIds[i] = rFemNeighbours[i]->Id();
            mFemIniNeighbourDelta[i] = initial_delta;
        }
    }//SetInitialFemContacts

    void SphericContinuumParticle::ContactAreaWeighting() //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbors of my neighbors.
    {
        double alpha = 1.0;
        //double external_sphere_area = 4 * Globals::Pi * GetRadius()*GetRadius();
        double effectiveVolumeRadius = EffectiveVolumeRadius();  //calculateEffectiveVolumeRadius
        double external_sphere_area = 4 * Globals::Pi * effectiveVolumeRadius * effectiveVolumeRadius;
        double total_equiv_area = 0.0;
        double total_mContIniNeighArea = 0.0;
        int cont_ini_neighbours_size = mContinuumInitialNeighborsSize;
        Vector& cont_ini_neigh_area = GetValue(NEIGHBOURS_CONTACT_AREAS);
        bool print_debug_files = false;

        for (int i = 0; i < cont_ini_neighbours_size; i++) {
            SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
            double other_radius = ini_cont_neighbour_iterator->GetRadius();
            double area = mContinuumConstitutiveLawArray[i]->CalculateContactArea(GetRadius(), other_radius, cont_ini_neigh_area); //This call fills the vector of areas only if the Constitutive Law wants.
            total_equiv_area += area;

        }
        if (print_debug_files) {
            std::ofstream outputfile("external_sphere_area-total_equiv_area.txt", std::ios_base::out | std::ios_base::app);
            outputfile << external_sphere_area << "  " << total_equiv_area << "\n";
            outputfile.close();
        }
        if (cont_ini_neighbours_size >= 6) {
            if (!IsSkin()) {
                AuxiliaryFunctions::CalculateAlphaFactor3D(cont_ini_neighbours_size, external_sphere_area, total_equiv_area, alpha);
                if (print_debug_files) {
                    std::ofstream outputfile("alpha.txt", std::ios_base::out | std::ios_base::app);
                    outputfile <<" "<<cont_ini_neighbours_size<<" "<<external_sphere_area<<" "<<total_equiv_area<<" "<< alpha << "\n";
                    outputfile.close();
                }
                for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                    cont_ini_neigh_area[i] = alpha * cont_ini_neigh_area[i];
                    total_mContIniNeighArea += cont_ini_neigh_area[i];
                } //for every neighbor
            }

            else {//skin sphere
                for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                    alpha = 1.00 * (1.40727)*(external_sphere_area / total_equiv_area)*((double(cont_ini_neighbours_size)) / 11.0);

                    cont_ini_neigh_area[i] = alpha * cont_ini_neigh_area[i];
                } //for every neighbor
            }
        }//if more than 6 neighbors
    }

    double SphericContinuumParticle::EffectiveVolumeRadius() {

        int cont_ini_neighbours_size = mContinuumInitialNeighborsSize;

        double effectiveVolumeRadiusSum = 0.0;

        for (int i = 0; i < cont_ini_neighbours_size; i++) {

            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            double other_radius = neighbour_iterator->GetRadius();
            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
            double distance = DEM_MODULUS_3(other_to_me_vect);
            effectiveVolumeRadiusSum += 0.5 * (distance + GetRadius() - other_radius);
        }

        double effectiveVolumeRadius = effectiveVolumeRadiusSum / cont_ini_neighbours_size;

        return effectiveVolumeRadius;
    }

    // TODO. Ignasi
    void SphericContinuumParticle::ComputeBallToBallInitialStiffness(SphericParticle::ParticleDataBuffer & data_buffer,
                                                            const ProcessInfo& r_process_info,
                                                            array_1d<double, 3>& r_nodal_stiffness_array,
                                                            array_1d<double, 3>& r_nodal_rotational_stiffness_array)
    {
        KRATOS_TRY

        NodeType& this_node = this->GetGeometry()[0];
        DEM_COPY_SECOND_TO_FIRST_3(data_buffer.mMyCoors, this_node)

        const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        Vector& cont_ini_neigh_area            = this->GetValue(NEIGHBOURS_CONTACT_AREAS);

        for (int i = 0; data_buffer.SetNextNeighbourOrExit(i); ++i) {

            if (mNeighbourElements[i] == NULL) continue;
            if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            data_buffer.mpOtherParticle = neighbour_iterator;

            noalias(data_buffer.mOtherToMeVector) = this->GetGeometry()[0].Coordinates() - data_buffer.mpOtherParticle->GetGeometry()[0].Coordinates();

            const double& other_radius = data_buffer.mpOtherParticle->GetRadius();

            data_buffer.mDistance = DEM_MODULUS_3(data_buffer.mOtherToMeVector);
            double radius_sum = GetRadius() + other_radius;

            double initial_delta = GetInitialDelta(i);

            double initial_dist = radius_sum - initial_delta;
            double indentation = initial_dist - data_buffer.mDistance;
            double myYoung = GetYoung();
            double myPoisson = GetPoisson();

            double kn_el = 0.0;
            double kt_el = 0.0;
            double DeltDisp[3] = {0.0};
            double RelVel[3] = {0.0};
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)

            // Getting neighbor properties
            double other_young = data_buffer.mpOtherParticle->GetYoung();
            double other_poisson = data_buffer.mpOtherParticle->GetPoisson();
            double equiv_poisson;
            if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
            else { equiv_poisson = 0.0; }

            double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
            double calculation_area = 0.0;
            const double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));

            if (i < (int)mContinuumInitialNeighborsSize) {
                mContinuumConstitutiveLawArray[i]->GetContactArea(GetRadius(), other_radius, cont_ini_neigh_area, i, calculation_area); //some Constitutive Laws get a value, some others calculate the value.
                mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area, this, neighbour_iterator, indentation);
            }

            EvaluateDeltaDisplacement(data_buffer, DeltDisp, RelVel, data_buffer.mLocalCoordSystem, data_buffer.mOldLocalCoordSystem, vel, delta_displ);

            double LocalStiffness[3]              = {0.0};
            double GlobalStiffness[3]             = {0.0};
            double LocalRotationalStiffness[3]       = {0.0};
            double GlobalRotationalStiffness[3]      = {0.0};

            if (i < (int)mContinuumInitialNeighborsSize) {

                LocalStiffness[0] = kt_el;
                LocalStiffness[1] = kt_el;
                LocalStiffness[2] = kn_el;

            } else if (indentation > 0.0) {

                mDiscontinuumConstitutiveLaw = pCloneDiscontinuumConstitutiveLawWithNeighbour(data_buffer.mpOtherParticle);
                mDiscontinuumConstitutiveLaw->InitializeContact(this,data_buffer.mpOtherParticle,indentation);

                const double Kt = mDiscontinuumConstitutiveLaw->mKt;
                const double Kn = mDiscontinuumConstitutiveLaw->mKn;

                LocalStiffness[0] = Kt;
                LocalStiffness[1] = Kt;
                LocalStiffness[2] = Kn;
            }

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                if (i < (int)mContinuumInitialNeighborsSize && mIniNeighbourFailureId[i] == 0) {

                    const double equivalent_radius = sqrt(calculation_area / Globals::Pi);
                    const double Inertia_I = 0.25 * Globals::Pi * equivalent_radius * equivalent_radius * equivalent_radius * equivalent_radius;
                    const double Inertia_J = 2.0 * Inertia_I; // This is the polar inertia
                    // Note. Ignasi: we do not reduce the rotational stiffness with the rotational_moment_coeff to avoid penalizing the delta time
                    // const double& rotational_moment_coeff = (*(mContinuumConstitutiveLawArray[i]->pGetProperties()))[ROTATIONAL_MOMENT_COEFFICIENT];
                    // LocalRotationalStiffness[0] = rotational_moment_coeff * equiv_young * Inertia_I / data_buffer.mDistance;
                    // LocalRotationalStiffness[1] = rotational_moment_coeff * equiv_young * Inertia_I / data_buffer.mDistance;
                    // LocalRotationalStiffness[2] = rotational_moment_coeff * equiv_young * Inertia_J / data_buffer.mDistance;
                    LocalRotationalStiffness[0] = equiv_young * Inertia_I / data_buffer.mDistance;
                    LocalRotationalStiffness[1] = equiv_young * Inertia_I / data_buffer.mDistance;
                    LocalRotationalStiffness[2] = equiv_young * Inertia_J / data_buffer.mDistance;
                }
            }

            GeometryFunctions::VectorLocal2Global(data_buffer.mLocalCoordSystem, LocalStiffness, GlobalStiffness);
            GeometryFunctions::VectorLocal2Global(data_buffer.mLocalCoordSystem, LocalRotationalStiffness, GlobalRotationalStiffness);

            // Make global stiffness positive before accumulation
            for(unsigned int j = 0; j<3; j++){
                if(GlobalStiffness[j] < 0.0) {
                    GlobalStiffness[j] = -GlobalStiffness[j];
                }
                if(GlobalRotationalStiffness[j] < 0.0) {
                    GlobalRotationalStiffness[j] = -GlobalRotationalStiffness[j];
                }
            }

            DEM_ADD_SECOND_TO_FIRST(r_nodal_stiffness_array, GlobalStiffness)
            DEM_ADD_SECOND_TO_FIRST(r_nodal_rotational_stiffness_array, GlobalRotationalStiffness)

        } // for each neighbor

        KRATOS_CATCH("")
    } //  ComputeBallToBallInitialStiffness

    void SphericContinuumParticle::ComputeBallToBallContactForce(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                const ProcessInfo& r_process_info,
                                                                array_1d<double, 3>& rElasticForce,
                                                                array_1d<double, 3>& rContactForce,
                                                                double& RollingResistance,
                                                                array_1d<double, 3>& r_nodal_stiffness_array,
                                                                array_1d<double, 3>& r_nodal_rotational_stiffness_array)
    {
        KRATOS_TRY

        NodeType& this_node = this->GetGeometry()[0];
        DEM_COPY_SECOND_TO_FIRST_3(data_buffer.mMyCoors, this_node)

        const int time_steps = r_process_info[TIME_STEPS];

        const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        Vector& cont_ini_neigh_area            = this->GetValue(NEIGHBOURS_CONTACT_AREAS);
        int NeighbourSize = mNeighbourElements.size();
        GetGeometry()[0].GetSolutionStepValue(NEIGHBOUR_SIZE) = NeighbourSize;

        for (int i = 0; data_buffer.SetNextNeighbourOrExit(i); ++i) {

            if (mNeighbourElements[i] == NULL) continue;
            if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            data_buffer.mpOtherParticle = neighbour_iterator;

            unsigned int neighbour_iterator_id = data_buffer.mpOtherParticle->Id();

            noalias(data_buffer.mOtherToMeVector) = this->GetGeometry()[0].Coordinates() - data_buffer.mpOtherParticle->GetGeometry()[0].Coordinates();

            const double& other_radius = data_buffer.mpOtherParticle->GetRadius();

            data_buffer.mDistance = DEM_MODULUS_3(data_buffer.mOtherToMeVector);
            double radius_sum = GetRadius() + other_radius;

            double initial_delta = GetInitialDelta(i);

            double initial_dist = radius_sum - initial_delta;
            double indentation = initial_dist - data_buffer.mDistance;
            double myYoung = GetYoung();
            double myPoisson = GetPoisson();

            double kn_el = 0.0;
            double kt_el = 0.0;
            double DeltDisp[3] = {0.0};
            double RelVel[3] = {0.0};
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)
            bool sliding = false;

            double contact_tau = 0.0;
            double contact_sigma = 0.0;
            double failure_criterion_state = 0.0;
            double acumulated_damage = 0.0;

            // Getting neighbor properties
            double other_young = data_buffer.mpOtherParticle->GetYoung();
            double other_poisson = data_buffer.mpOtherParticle->GetPoisson();
            double equiv_poisson;
            if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
            else { equiv_poisson = 0.0; }

            double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
            double calculation_area = 0.0;
            const double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));

            if (i < (int)mContinuumInitialNeighborsSize) {
                mContinuumConstitutiveLawArray[i]->GetContactArea(GetRadius(), other_radius, cont_ini_neigh_area, i, calculation_area); //some Constitutive Laws get a value, some others calculate the value.
                mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area, this, neighbour_iterator, indentation);
            }

            EvaluateDeltaDisplacement(data_buffer, DeltDisp, RelVel, data_buffer.mLocalCoordSystem, data_buffer.mOldLocalCoordSystem, vel, delta_displ);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                RelativeDisplacementAndVelocityOfContactPointDueToRotationQuaternion(DeltDisp, RelVel, data_buffer.mOldLocalCoordSystem, other_radius, data_buffer.mDt, ang_vel, neighbour_iterator);
            }

            RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(r_process_info, DeltDisp, RelVel, data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, neighbour_iterator);

            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3] = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double LocalElasticExtraContactForce[3] = {0.0};
            double GlobalElasticContactForce[3] = {0.0};
            double GlobalElasticExtraContactForce[3] = {0.0};
            double TotalGlobalElasticContactForce[3] = {0.0};
            double OldLocalElasticContactForce[3] = {0.0};
            double LocalStiffness[3]              = {0.0};
            double GlobalStiffness[3]             = {0.0};
            double LocalRotationalStiffness[3]    = {0.0};
            double GlobalRotationalStiffness[3]   = {0.0};

            //FilterNonSignificantDisplacements(DeltDisp, RelVel, indentation);

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltDisp, LocalDeltDisp);

            RotateOldContactForces(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, mNeighbourElasticContactForces[i]);
            RotateOldContactForces(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, mNeighbourElasticExtraContactForces[i]);

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce);

            GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i][0];
            GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i][1];
            GlobalElasticContactForce[2] = mNeighbourElasticContactForces[i][2];

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); //TODO: can we remove this? We should overwrite LocalElasticContactForce afterwards

            double ViscoDampingLocalContactForce[3] = {0.0};
            double equiv_visco_damp_coeff_normal;
            double equiv_visco_damp_coeff_tangential;
            double ElasticLocalRotationalMoment[3] = {0.0};
            double ViscoLocalRotationalMoment[3] = {0.0};
            double cohesive_force =  0.0;
            double LocalRelVel[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, RelVel, LocalRelVel);

            if (i < (int)mContinuumInitialNeighborsSize) {

                mContinuumConstitutiveLawArray[i]->CheckFailure(i, this, neighbour_iterator);

                mContinuumConstitutiveLawArray[i]->CalculateForces(r_process_info,
                                                                OldLocalElasticContactForce,
                                                                LocalElasticContactForce,
                                                                LocalElasticExtraContactForce,
                                                                data_buffer.mLocalCoordSystem,
                                                                LocalDeltDisp,
                                                                kn_el,
                                                                kt_el,
                                                                contact_sigma,
                                                                contact_tau,
                                                                failure_criterion_state,
                                                                equiv_young,
                                                                equiv_shear,
                                                                indentation,
                                                                calculation_area,
                                                                acumulated_damage,
                                                                this,
                                                                neighbour_iterator,
                                                                i,
                                                                r_process_info[TIME_STEPS],
                                                                sliding,
                                                                equiv_visco_damp_coeff_normal,
                                                                equiv_visco_damp_coeff_tangential,
                                                                LocalRelVel,
                                                                ViscoDampingLocalContactForce);

                // Estimate updated local stiffness
                if (std::abs(LocalDeltDisp[0]) > std::numeric_limits<double>::epsilon()){
                    LocalStiffness[0] = -(LocalElasticContactForce[0]-OldLocalElasticContactForce[0]) / LocalDeltDisp[0];
                } else{
                    LocalStiffness[0] = kt_el;
                }
                if (std::abs(LocalDeltDisp[1]) > std::numeric_limits<double>::epsilon()){
                    LocalStiffness[1] = -(LocalElasticContactForce[1]-OldLocalElasticContactForce[1]) / LocalDeltDisp[1];
                } else{
                    LocalStiffness[1] = kt_el;
                }
                if (std::abs(indentation) > std::numeric_limits<double>::epsilon()){
                    LocalStiffness[2] = LocalElasticContactForce[2] / indentation;
                } else{
                    LocalStiffness[2] = kn_el;
                }

            } else if (indentation > 0.0) {
                const double previous_indentation = indentation + LocalDeltDisp[2];
                mDiscontinuumConstitutiveLaw = pCloneDiscontinuumConstitutiveLawWithNeighbour(data_buffer.mpOtherParticle);
                mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce,
                        LocalDeltDisp, LocalRelVel, indentation, previous_indentation,
                        ViscoDampingLocalContactForce, cohesive_force, this, data_buffer.mpOtherParticle, sliding, data_buffer.mLocalCoordSystem);
                
                // Estimate updated local stiffness
                mDiscontinuumConstitutiveLaw->InitializeContact(this,data_buffer.mpOtherParticle,indentation);
                const double Kt = mDiscontinuumConstitutiveLaw->mKt;
                const double Kn = mDiscontinuumConstitutiveLaw->mKn;
                if (std::abs(LocalDeltDisp[0]) > std::numeric_limits<double>::epsilon()){
                    LocalStiffness[0] = -(LocalElasticContactForce[0]-OldLocalElasticContactForce[0]) / LocalDeltDisp[0];
                } else{
                    LocalStiffness[0] = Kt;
                }
                if (std::abs(LocalDeltDisp[1]) > std::numeric_limits<double>::epsilon()){
                    LocalStiffness[1] = -(LocalElasticContactForce[1]-OldLocalElasticContactForce[1]) / LocalDeltDisp[1];
                } else{
                    LocalStiffness[1] = Kt;
                }
                if (std::abs(indentation) > std::numeric_limits<double>::epsilon()){
                    LocalStiffness[2] = LocalElasticContactForce[2] / indentation;
                } else{
                    LocalStiffness[2] = Kn;
                }
            } else { //Not bonded and no idata_buffer.mpOtherParticlendentation
                LocalElasticContactForce[0] = 0.0;      LocalElasticContactForce[1] = 0.0;      LocalElasticContactForce[2] = 0.0;
                ViscoDampingLocalContactForce[0] = 0.0; ViscoDampingLocalContactForce[1] = 0.0; ViscoDampingLocalContactForce[2] = 0.0;
                cohesive_force= 0.0;
            }

            // Transforming to global forces and adding up
            double LocalContactForce[3] = {0.0};
            double GlobalContactForce[3] = {0.0};

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < (int)mContinuumInitialNeighborsSize)) { // We leave apart the discontinuum neighbors (the same for the walls). The neighbor would not be able to do the same if we activate it.
                mContinuumConstitutiveLawArray[i]->AddPoissonContribution(equiv_poisson, data_buffer.mLocalCoordSystem, LocalElasticContactForce[2], calculation_area, mSymmStressTensor, this, neighbour_iterator, r_process_info, i, indentation);
            }

            array_1d<double, 3> other_ball_to_ball_forces(3,0.0);
            ComputeOtherBallToBallForces(other_ball_to_ball_forces);

            AddUpForcesAndProject(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, LocalContactForce, LocalElasticContactForce, LocalElasticExtraContactForce, GlobalContactForce,
                                  GlobalElasticContactForce, GlobalElasticExtraContactForce, TotalGlobalElasticContactForce, ViscoDampingLocalContactForce, 0.0, other_ball_to_ball_forces, rElasticForce, rContactForce, i, r_process_info); //TODO: replace the 0.0 with an actual cohesive force for discontinuum neighbours

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                ComputeMoments(LocalContactForce[2], TotalGlobalElasticContactForce, RollingResistance, data_buffer.mLocalCoordSystem[2], data_buffer.mpOtherParticle, indentation, i);
                if (i < (int)mContinuumInitialNeighborsSize && mIniNeighbourFailureId[i] == 0) {
                    mContinuumConstitutiveLawArray[i]->ComputeParticleRotationalMoments(this, neighbour_iterator, equiv_young, data_buffer.mDistance, calculation_area,
                                                                                        data_buffer.mLocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment, equiv_poisson, indentation);

                    // Estimate updated local rotational stiffness
                    double LocalDeltaRotatedAngle[3]    = {0.0};
                    array_1d<double, 3> GlobalDeltaRotatedAngle;
                    noalias(GlobalDeltaRotatedAngle) = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE) - neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
                    GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, GlobalDeltaRotatedAngle, LocalDeltaRotatedAngle);
                    double aux = (this->GetRadius() + neighbour_iterator->GetRadius()) / data_buffer.mDistance;
                    array_1d<double, 3> LocalEffDeltaRotatedAngle;
                    LocalEffDeltaRotatedAngle[0] = LocalDeltaRotatedAngle[0] * aux;
                    LocalEffDeltaRotatedAngle[1] = LocalDeltaRotatedAngle[1] * aux;
                    LocalEffDeltaRotatedAngle[2] = LocalDeltaRotatedAngle[2] * aux;
                    const double equivalent_radius = sqrt(calculation_area / Globals::Pi);
                    const double Inertia_I = 0.25 * Globals::Pi * equivalent_radius * equivalent_radius * equivalent_radius * equivalent_radius;
                    const double Inertia_J = 2.0 * Inertia_I; // This is the polar inertia
                    const double& rotational_moment_coeff = (*(mContinuumConstitutiveLawArray[i]->pGetProperties()))[ROTATIONAL_MOMENT_COEFFICIENT];

                    if (std::abs(LocalEffDeltaRotatedAngle[0]) > std::numeric_limits<double>::epsilon()){
                        // LocalRotationalStiffness[0] = -ElasticLocalRotationalMoment[0] / LocalEffDeltaRotatedAngle[0];
                        LocalRotationalStiffness[0] = -ElasticLocalRotationalMoment[0] / (LocalEffDeltaRotatedAngle[0] * rotational_moment_coeff);
                    } else{
                        // LocalRotationalStiffness[0] = rotational_moment_coeff * equiv_young * Inertia_I / data_buffer.mDistance;
                        LocalRotationalStiffness[0] = equiv_young * Inertia_I / data_buffer.mDistance;
                    }
                    if (std::abs(LocalEffDeltaRotatedAngle[1]) > std::numeric_limits<double>::epsilon()){
                        // LocalRotationalStiffness[1] = -ElasticLocalRotationalMoment[1] / LocalEffDeltaRotatedAngle[1];
                        LocalRotationalStiffness[1] = -ElasticLocalRotationalMoment[1] / (LocalEffDeltaRotatedAngle[1] * rotational_moment_coeff);
                    } else{
                        // LocalRotationalStiffness[1] = rotational_moment_coeff * equiv_young * Inertia_I / data_buffer.mDistance;
                        LocalRotationalStiffness[1] = equiv_young * Inertia_I / data_buffer.mDistance;
                    }
                    if (std::abs(LocalEffDeltaRotatedAngle[2]) > std::numeric_limits<double>::epsilon()){
                        // LocalRotationalStiffness[2] = -ElasticLocalRotationalMoment[2] / LocalEffDeltaRotatedAngle[2];
                        LocalRotationalStiffness[2] = -ElasticLocalRotationalMoment[2] / (LocalEffDeltaRotatedAngle[2] * rotational_moment_coeff);
                    } else{
                        // LocalRotationalStiffness[2] = rotational_moment_coeff * equiv_young * Inertia_J / data_buffer.mDistance;
                        LocalRotationalStiffness[2] = equiv_young * Inertia_J / data_buffer.mDistance;
                    }
                }

                AddUpMomentsAndProject(data_buffer.mLocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
            }

            // Accumulate GlobalStiffness
            GeometryFunctions::VectorLocal2Global(data_buffer.mLocalCoordSystem, LocalStiffness, GlobalStiffness);
            GeometryFunctions::VectorLocal2Global(data_buffer.mLocalCoordSystem, LocalRotationalStiffness, GlobalRotationalStiffness);
            // Make global stiffness positive before accumulation
            for(unsigned int j = 0; j<3; j++){
                if(GlobalStiffness[j] < 0.0) {
                    GlobalStiffness[j] = -GlobalStiffness[j];
                }
                if(GlobalRotationalStiffness[j] < 0.0) {
                    GlobalRotationalStiffness[j] = -GlobalRotationalStiffness[j];
                }
            }
            DEM_ADD_SECOND_TO_FIRST(r_nodal_stiffness_array, GlobalStiffness)
            DEM_ADD_SECOND_TO_FIRST(r_nodal_rotational_stiffness_array, GlobalRotationalStiffness)

            if (r_process_info[CONTACT_MESH_OPTION] == 1 && (i < (int)mContinuumInitialNeighborsSize) && this->Id() < neighbour_iterator_id) {
                double total_local_elastic_contact_force[3] = {0.0};
                total_local_elastic_contact_force[0] = LocalElasticContactForce[0] + LocalElasticExtraContactForce[0];
                total_local_elastic_contact_force[1] = LocalElasticContactForce[1] + LocalElasticExtraContactForce[1];
                total_local_elastic_contact_force[2] = LocalElasticContactForce[2] + LocalElasticExtraContactForce[2];
                CalculateOnContinuumContactElements(i, total_local_elastic_contact_force, contact_sigma, contact_tau, failure_criterion_state, acumulated_damage, time_steps);
            }

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) /*&& (i < mContinuumInitialNeighborsSize)*/) {
                AddNeighbourContributionToStressTensor(r_process_info, TotalGlobalElasticContactForce, data_buffer.mLocalCoordSystem[2], data_buffer.mDistance, radius_sum, this);
            }

            AddContributionToRepresentativeVolume(data_buffer.mDistance, radius_sum, calculation_area);

            ComputeForceWithNeighbourFinalOperations();

            /*if (i < mContinuumInitialNeighborsSize) {
                DEM_COPY_SECOND_TO_FIRST_3(mArrayOfDeltaDisplacements[i], DeltDisp);
            }*/
        } // for each neighbor

        ComputeBrokenBondsRatio();

        KRATOS_CATCH("")
    } //  ComputeBallToBallContactForce

    void SphericContinuumParticle::ComputeRollingResistance(double& RollingResistance, const double& NormalLocalContactForce, const double& equiv_rolling_friction_coeff, const unsigned int i) {

        KRATOS_TRY

        if ((i >= mContinuumInitialNeighborsSize) || mIniNeighbourFailureId[i]) {
            SphericParticle::ComputeRollingResistance(RollingResistance, NormalLocalContactForce, equiv_rolling_friction_coeff, i);
        }

        KRATOS_CATCH("")
    } // ComputeRollingResistance

    void SphericContinuumParticle::ComputeForceWithNeighbourFinalOperations(){}

    void SphericContinuumParticle::ComputeBrokenBondsRatio() {

        int BrokenBondsCounter = 0.0;

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {

            if (mNeighbourElements[i] == NULL) continue;
            //if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;

            //SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            //unsigned int neighbour_id = neighbour_iterator->Id();
            if ((i < mContinuumInitialNeighborsSize) && mIniNeighbourFailureId[i] > 0) {
            //if ((i < mContinuumInitialNeighborsSize) && this->Id() < neighbour_id && mIniNeighbourFailureId[i] > 0) {
                ++BrokenBondsCounter;
            }

            //int NeighbourSize = mNeighbourElements.size();
            //GetGeometry()[0].GetSolutionStepValue(NEIGHBOUR_SIZE) = NeighbourSize;
            double NeighbourRatio = 0.0;

            if (mContinuumInitialNeighborsSize) NeighbourRatio = BrokenBondsCounter / mContinuumInitialNeighborsSize;

            GetGeometry()[0].GetSolutionStepValue(NEIGHBOUR_RATIO) = NeighbourRatio;
        }
    }

    void SphericContinuumParticle::FilterNonSignificantDisplacements(double DeltDisp[3], //IN GLOBAL AXES
                                                                     double RelVel[3], //IN GLOBAL AXES
                                                                     double& indentation) {

        if (fabs(DeltDisp[0]) < 1e-15){DeltDisp[0] = 0.0;}
        if (fabs(DeltDisp[1]) < 1e-15){DeltDisp[1] = 0.0;}
        if (fabs(DeltDisp[2]) < 1e-15){DeltDisp[2] = 0.0;}
        if (fabs(RelVel[0]) < 1e-15){RelVel[0] = 0.0;}
        if (fabs(RelVel[1]) < 1e-15){RelVel[1] = 0.0;}
        if (fabs(RelVel[2]) < 1e-15){RelVel[2] = 0.0;}
        if (fabs(indentation) < 1e-15){indentation = 0.0;}
    }

    void SphericContinuumParticle::AddContributionToRepresentativeVolume(const double distance,
                                                                         const double radius_sum,
                                                                         const double contact_area) {
        KRATOS_TRY

        double gap = distance - radius_sum;
        double real_distance = GetInteractionRadius() + 0.5 * gap;
        //double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
        //rRepresentative_Volume += 0.33333333333333 * (real_distance * contact_area);
        mPartialRepresentativeVolume += 0.33333333333333 * real_distance * contact_area;

        KRATOS_CATCH("")
    }

    /*
    void SphericContinuumParticle::CorrectRepresentativeVolume(double& rRepresentative_Volume, bool& is_smaller_than_sphere) {

        KRATOS_TRY

        SphericParticle::CorrectRepresentativeVolume(rRepresentative_Volume, is_smaller_than_sphere);

        //if (*mSkinSphere && is_smaller_than_sphere) rRepresentative_Volume *= 1.5; // This is the quotient between the volume of the cylinder circumscribed about a sphere and the latter
                                                                                   // So the minimum volume for a skin sphere is that of a cylinder
        //if (*mSkinSphere) rRepresentative_Volume *= 0.5;

        KRATOS_CATCH("")
    }
    */

    void SphericContinuumParticle::FinalizeSolutionStep(const ProcessInfo& r_process_info) {
        KRATOS_TRY

        SphericParticle::FinalizeSolutionStep(r_process_info);

        //Update sphere mass and inertia taking into account the real volume of the represented volume:
        SetMass(this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) * GetDensity());
        if (this->Is(DEMFlags::HAS_ROTATION) ){
            GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = CalculateMomentOfInertia();
        }

        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::GetStressTensorFromNeighbourStep1(){

        Set(DEMFlags::COPIED_STRESS_TENSOR, false);
        Set(DEMFlags::COPIED_STRESS_TENSOR2,false);
        if(!IsSkin()) return;

        for (unsigned int i=0; i<mNeighbourElements.size(); i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericContinuumParticle* p_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            if(!p_neighbour->IsSkin()) {
                *(mStressTensor) = *(p_neighbour->mStressTensor);
                *(mSymmStressTensor) = *(p_neighbour->mSymmStressTensor);
                Set(DEMFlags::COPIED_STRESS_TENSOR,true);
                break;
            }
        }
    }

    void SphericContinuumParticle::GetStressTensorFromNeighbourStep2(){

        if(!IsSkin()) return;
        if(Is(DEMFlags::COPIED_STRESS_TENSOR)) return;

        for (unsigned int i=0; i<mNeighbourElements.size(); i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericContinuumParticle* p_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            if(p_neighbour->Is(DEMFlags::COPIED_STRESS_TENSOR)) {
                *(mStressTensor) = *(p_neighbour->mStressTensor);
                *(mSymmStressTensor) = *(p_neighbour->mSymmStressTensor);
                Set(DEMFlags::COPIED_STRESS_TENSOR2,true);
                break;
            }
        }
    }

    void SphericContinuumParticle::GetStressTensorFromNeighbourStep3(){

        if(!IsSkin()) return;
        if(Is(DEMFlags::COPIED_STRESS_TENSOR)) return;
        if(Is(DEMFlags::COPIED_STRESS_TENSOR2)) return;

        for (unsigned int i=0; i<mNeighbourElements.size(); i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericContinuumParticle* p_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            if(p_neighbour->Is(DEMFlags::COPIED_STRESS_TENSOR2)) {
                *(mStressTensor) = *(p_neighbour->mStressTensor);
                *(mSymmStressTensor) = *(p_neighbour->mSymmStressTensor);
                break;
            }
        }
    }

    void SphericContinuumParticle::MarkNewSkinParticlesDueToBreakage() {

        KRATOS_TRY

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {

            if (mNeighbourElements[i] == NULL || mIniNeighbourFailureId[i]) {
                *mSkinSphere = 1.0;
                break;
            }
        }

        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::ReorderAndRecoverInitialPositionsAndFilter(std::vector<SphericParticle*>& temp_neighbour_elements) {

        KRATOS_TRY

        unsigned int current_neighbors_size = mNeighbourElements.size();
        unsigned int initial_neighbors_size = mIniNeighbourIds.size();
        temp_neighbour_elements.resize(initial_neighbors_size);

        for (unsigned int i = 0; i < initial_neighbors_size; i++) {
            temp_neighbour_elements[i] = NULL;
        }

        // Loop over current neighbors
        for (unsigned int i = 0; i < current_neighbors_size; i++) {
            SphericParticle* i_neighbour = mNeighbourElements[i];
            bool found = false;
            // Loop over initial neighbors
            for (unsigned int kk = 0; kk < initial_neighbors_size; kk++) {

                if (static_cast<int>(i_neighbour->Id()) == mIniNeighbourIds[kk]) {
                    temp_neighbour_elements[kk] = i_neighbour;
                    found = true;
                    break;
                }
            }

            if (!found) {
                double other_radius = i_neighbour->GetInteractionRadius();
                double radius_sum = GetInteractionRadius() + other_radius;
                array_1d<double, 3> other_to_me_vect;
                noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - i_neighbour->GetGeometry()[0].Coordinates();
                double distance = DEM_MODULUS_3(other_to_me_vect);
                double indentation = radius_sum - distance;

                if (indentation > 0.0) {
                    temp_neighbour_elements.push_back(i_neighbour);
                }
            }
        }

        mNeighbourElements.swap(temp_neighbour_elements);

        if (mBondElements.size()) {
            for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
                if (mNeighbourElements[i] == NULL) {
                    mBondElements[i] = NULL;
                }
                if ((mNeighbourElements[i] == NULL) && (mIniNeighbourFailureId[i] == 0)) {
                    mIniNeighbourFailureId[i] = 6; // Breakage due to search (neighbor not broken was not found during search)
                }
            }
        }

        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::UpdateContinuumNeighboursVector(const ProcessInfo& r_process_info) {}

    void SphericContinuumParticle::ReorderFEMneighbours() {

        KRATOS_TRY

        unsigned int current_neighbors_size = mNeighbourRigidFaces.size();
        unsigned int initial_neighbors_size = mFemIniNeighbourIds.size();

        std::vector<DEMWall*> temp_neighbour_elements(initial_neighbors_size, nullptr);
        std::vector<array_1d<double, 4> > temp_neighbours_weights(initial_neighbors_size, ZeroVector(4));
        std::vector<int> temp_neighbours_contact_types(initial_neighbors_size, 0);

        // Loop over current neighbors
        for (unsigned int i = 0; i < current_neighbors_size; i++) {
            DEMWall* i_neighbour = mNeighbourRigidFaces[i];
            bool found = false;
            // Loop over initial neighbors
            for (unsigned int k = 0; k < initial_neighbors_size; k++) {
                if (static_cast<int>(i_neighbour->Id()) == mFemIniNeighbourIds[k]) {
                    temp_neighbour_elements[k] = i_neighbour;
                    temp_neighbours_weights[k] = mContactConditionWeights[i];
                    temp_neighbours_contact_types[k] = mContactConditionContactTypes[i];
                    found = true;
                    break;
                }
            }

            if (!found) {
                temp_neighbour_elements.push_back(i_neighbour);
                temp_neighbours_weights.push_back(mContactConditionWeights[i]);
                temp_neighbours_contact_types.push_back(mContactConditionContactTypes[i]);
            }
        }

        mNeighbourRigidFaces.swap(temp_neighbour_elements);
        mContactConditionWeights.swap(temp_neighbours_weights);
        mContactConditionContactTypes.swap(temp_neighbours_contact_types);

        KRATOS_CATCH("")
    }

    double SphericContinuumParticle::CalculateMaxSearchDistance(const bool has_mpi, const ProcessInfo& r_process_info) {

        KRATOS_TRY

        double max_local_search = 0.0;

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericContinuumParticle* r_continuum_ini_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            double search_dist = mContinuumConstitutiveLawArray[i]->LocalMaxSearchDistance(i, this, r_continuum_ini_neighbour);
            if (search_dist > max_local_search) max_local_search = search_dist;
        }

        return max_local_search;

        KRATOS_CATCH("")
    }

    bool SphericContinuumParticle::OverlappedParticleRemoval() {

        KRATOS_TRY

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
            double other_radius = ini_cont_neighbour_iterator->GetRadius();

            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - mNeighbourElements[i]->GetGeometry()[0].Coordinates();
            double distance = DEM_MODULUS_3(other_to_me_vect);

            double alpha = 1.0; // alpha = 1.0 means that the particle is completely inside another

            if (GetRadius() + distance < alpha * other_radius) {
                this->Set(TO_ERASE, true);
                return true;
            }
        }
        return false;

        KRATOS_CATCH("")
    }

    double SphericContinuumParticle::CalculateLocalMaxPeriod(const bool has_mpi, const ProcessInfo& r_process_info) {

        KRATOS_TRY

        double max_sqr_period = 0.0;

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            SphericContinuumParticle* r_continuum_ini_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            double sqr_period = mContinuumConstitutiveLawArray[i]->LocalPeriod(i, this, r_continuum_ini_neighbour);
            if (sqr_period > max_sqr_period) { (max_sqr_period = sqr_period); }
        }
        for (unsigned int i = mContinuumInitialNeighborsSize; i < mNeighbourElements.size(); i++) {
            mDiscontinuumConstitutiveLaw = pCloneDiscontinuumConstitutiveLawWithNeighbour(mNeighbourElements[i]);
            double sqr_period_discontinuum = mDiscontinuumConstitutiveLaw->LocalPeriod(i, this, mNeighbourElements[i]);
            if (sqr_period_discontinuum > max_sqr_period) { (max_sqr_period = sqr_period_discontinuum); }
        }
        return max_sqr_period;

        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::CalculateMeanContactArea(const bool has_mpi, const ProcessInfo& r_process_info) {

        KRATOS_TRY

        Vector& cont_ini_neigh_area = GetValue(NEIGHBOURS_CONTACT_AREAS);

        if (!cont_ini_neigh_area.size()) return; //TODO: ugly fix //This means that for this case the areas are not being saved (because the constitutive law is not filling this vector)

        for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {

            SphericContinuumParticle* r_continuum_ini_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            if (r_continuum_ini_neighbour == NULL) continue; //The initial neighbor was deleted at some point in time!!

            if (this->Id() > r_continuum_ini_neighbour->Id()) continue; // The Sphere with a lower Id will do the job only

            Vector& NeighbourContIniNeighArea = r_continuum_ini_neighbour->GetValue(NEIGHBOURS_CONTACT_AREAS);

            int index_of_the_neighbour_that_is_me = -1;

            for (unsigned int j = 0; j <  NeighbourContIniNeighArea.size(); j++) {
                DenseVector<int>& vector_of_ids_of_neighbours_of_this_neighbour = r_continuum_ini_neighbour->GetValue(NEIGHBOUR_IDS);
                if ((int) this->Id() == vector_of_ids_of_neighbours_of_this_neighbour[j]) index_of_the_neighbour_that_is_me = (int)j;
            }

            if (index_of_the_neighbour_that_is_me == -1) {
                std::string message = "An element (Id " + std::to_string(this->Id()) + ") found a neighbor (had contact area) but the neighbor (Id " \
                                                        + std::to_string(r_continuum_ini_neighbour->Id()) + ") did not have area for that element  ";

                KRATOS_ERROR << message << std::endl;
            }

            bool neigh_is_skin = r_continuum_ini_neighbour->IsSkin();
            if ((IsSkin() && neigh_is_skin) || (!IsSkin() && !neigh_is_skin)) { //both skin or both inner
                double mean_area =  0.5 * (cont_ini_neigh_area[i] + NeighbourContIniNeighArea[index_of_the_neighbour_that_is_me]);
                cont_ini_neigh_area[i] = mean_area;
                NeighbourContIniNeighArea[index_of_the_neighbour_that_is_me] = mean_area;
            }
            else if (!IsSkin() && neigh_is_skin) {//we will store both the same only coming from the inner to the skin.
                NeighbourContIniNeighArea[index_of_the_neighbour_that_is_me] = cont_ini_neigh_area[i];
            }
            else {
                cont_ini_neigh_area[i] = NeighbourContIniNeighArea[index_of_the_neighbour_that_is_me];
            }
        } //loop neighbors

        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {

        KRATOS_TRY

        if (rVariable == DELTA_TIME) {
            double coeff = r_process_info[NODAL_MASS_COEFF];
            double mass = GetMass();
            KRATOS_ERROR_IF(coeff > 1.0) << "The coefficient assigned for virtual mass is larger than one, virtual_mass_coeff= "<< coeff << std::endl;
            else if ((coeff == 1.0) && (r_process_info[VIRTUAL_MASS_OPTION])) { Output = 9.0E09; }
            else {
                if (r_process_info[VIRTUAL_MASS_OPTION]) { mass /= 1 - coeff; }
                double K = GetYoung() * Globals::Pi * GetRadius();
                Output = 0.34 * sqrt(mass / K);
                if (r_process_info[ROTATION_OPTION] == 1) { Output = Output * 0.5; } //factor for critical time step when rotation is allowed.
            }
            return;
        }//CRITICAL DELTA CALCULATION

        SphericParticle::Calculate(rVariable, Output, r_process_info);

        KRATOS_CATCH("")
    }//Calculate

    void SphericContinuumParticle::Initialize(const ProcessInfo& r_process_info){

        SphericParticle::Initialize(r_process_info);

        SetValue(NEIGHBOURS_CONTACT_AREAS, Vector());

        mSkinSphere     = &(this->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE));
        mContinuumGroup = this->GetGeometry()[0].FastGetSolutionStepValue(COHESIVE_GROUP);
    }

    double SphericContinuumParticle::GetInitialDelta(int index) {
        if (index < (int) mIniNeighbourDelta.size()) return mIniNeighbourDelta[index];
        else  return 0.0;
    }

    double SphericContinuumParticle::GetInitialDeltaWithFEM(int index) {
        if(index< (int) mFemIniNeighbourDelta.size()) return mFemIniNeighbourDelta[index];
        else return 0.0;
    }

    void SphericContinuumParticle::CalculateOnContinuumContactElements(size_t i, double LocalElasticContactForce[3], double contact_sigma, double contact_tau, double failure_criterion_state, double acumulated_damage, int time_steps) {

        KRATOS_TRY
        if (!mBondElements.size()) return; // we skip this function if the vector of bonds hasn't been filled yet.
        ParticleContactElement* bond = mBondElements[i];
        if (bond == NULL) return; //This bond was never created (happens in some MPI cases, see CreateContactElements() in explicit_solve_continumm.h)

        bond->mLocalContactForce[0] = LocalElasticContactForce[0];
        bond->mLocalContactForce[1] = LocalElasticContactForce[1];
        bond->mLocalContactForce[2] = LocalElasticContactForce[2];
        bond->mContactSigma = contact_sigma;
        bond->mContactTau = contact_tau;
        bond->mContactFailure = mIniNeighbourFailureId[i];
        bond->mFailureCriterionState = failure_criterion_state;

        if ((time_steps == 0) || (acumulated_damage > bond->mUnidimendionalDamage)) {
            bond->mUnidimendionalDamage = acumulated_damage;
        }
        KRATOS_CATCH("")
    }

} // namespace Kratos

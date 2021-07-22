
proc WriteMaterialsDEM { basename dir problemtypedir } {

    set filename [file join $dir MaterialsDEM.json]
    set FileVar  [open $filename w]

    puts $FileVar "\{"
    puts $FileVar "    \"materials\" : \[ \{ "
    puts $FileVar "        \"material_name\": \"mat1\","
    puts $FileVar "        \"material_id\": 1,"
    puts $FileVar "        \"Variables\": \{"
    puts $FileVar "            \"PARTICLE_DENSITY\": [GiD_AccessValue get gendata Density],"
    puts $FileVar "            \"YOUNG_MODULUS\": [GiD_AccessValue get gendata Young_Modulus],"
    puts $FileVar "            \"POISSON_RATIO\": [GiD_AccessValue get gendata Poisson_Ratio]"
    puts $FileVar "        \}"
    puts $FileVar "    \},\{ "
    puts $FileVar "        \"material_name\": \"mat2\","
    puts $FileVar "        \"material_id\": 2,"
    puts $FileVar "        \"Variables\": \{"
    puts $FileVar "            \"YOUNG_MODULUS\": [GiD_AccessValue get gendata Young_Modulus],"
    puts $FileVar "            \"POISSON_RATIO\": [GiD_AccessValue get gendata Poisson_Ratio],"
    puts $FileVar "            \"COMPUTE_WEAR\": false"
    puts $FileVar "        \}"
    puts $FileVar "    \}\],"
    puts $FileVar "    \"material_relations\" : \[ \{"
    puts $FileVar "        \"material_names_list\" : \[ \"mat1\", \"mat1\" \],"
    puts $FileVar "        \"material_ids_list\": \[1, 1\],"
    puts $FileVar "        \"Variables\":\{"
    puts $FileVar "            \"COEFFICIENT_OF_RESTITUTION\": [GiD_AccessValue get gendata Coefficion_of_Restitution],"
    puts $FileVar "            \"STATIC_FRICTION\": [GiD_AccessValue get gendata Static_Friction],"
    puts $FileVar "            \"DYNAMIC_FRICTION\":[GiD_AccessValue get gendata Dynamic_Friction],"
    puts $FileVar "            \"FRICTION_DECAY\": 500,"
    puts $FileVar "            \"ROLLING_FRICTION\": [GiD_AccessValue get gendata Rolling_Friction],"
    puts $FileVar "            \"ROLLING_FRICTION_WITH_WALLS\": [GiD_AccessValue get gendata Rolling_Friction],"
    puts $FileVar "            \"DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME\": \"DEM_D_Linear_viscous_Coulomb\""
    puts $FileVar "        \}"
    puts $FileVar "    \},\{"
    puts $FileVar "        \"material_names_list\":\[\"mat1\", \"mat2\"\],"
    puts $FileVar "        \"material_ids_list\":\[1, 2\],"
    puts $FileVar "        \"Variables\":\{"
    puts $FileVar "            \"COEFFICIENT_OF_RESTITUTION\": 0.2,"
    puts $FileVar "            \"STATIC_FRICTION\": [GiD_AccessValue get gendata Static_Friction],"
    puts $FileVar "            \"DYNAMIC_FRICTION\": [GiD_AccessValue get gendata Dynamic_Friction],"
    puts $FileVar "            \"FRICTION_DECAY\": 500,"
    puts $FileVar "            \"ROLLING_FRICTION\": [GiD_AccessValue get gendata Rolling_Friction],"
    puts $FileVar "            \"ROLLING_FRICTION_WITH_WALLS\": 0.01,"
    puts $FileVar "            \"SEVERITY_OF_WEAR\": 0.001,"
    puts $FileVar "            \"IMPACT_WEAR_SEVERITY\": 0.001,"
    puts $FileVar "            \"BRINELL_HARDNESS\": 200.0,"
    puts $FileVar "            \"DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME\": \"DEM_D_Linear_viscous_Coulomb\""
    puts $FileVar "        \}"
    puts $FileVar "    \}\],"
    puts $FileVar "    \"material_assignation_table\" :\["
    puts $FileVar "    \]"
    puts $FileVar "\}"

    close $FileVar
}
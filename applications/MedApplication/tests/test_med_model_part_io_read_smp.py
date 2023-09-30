import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.MedApplication as KratosMed

from itertools import chain
from pathlib import Path
from typing import Any, Dict, List, Set


def GetMedPath(med_path, med_name="mesh.med"):
    return Path(__file__).absolute().parent / "med_files" / med_path / med_name


class TestMedModelPartIOReadSubModelPart(KratosUnittest.TestCase):

    def test_cube_with_groups(self):
        model: KM.Model = KM.Model()
        model_part: KM.ModelPart = model.CreateModelPart("test")
        KratosMed.MedModelPartIO(GetMedPath("cube_with_groups", "cube.med")).ReadModelPart(model_part)

        self._basic_checks(model_part)

        self.assertEqual(model_part.NumberOfNodes(), 130)
        self.assertEqual(model_part.NumberOfGeometries(), 644)

        self.assertEqual(model_part.NumberOfSubModelParts(), 2)
        self.assertTrue(model_part.HasSubModelPart("interface"))
        self.assertTrue(model_part.HasSubModelPart("interface_nodes"))

        smp_interface: KM.ModelPart = model_part.GetSubModelPart("interface")
        smp_interface_nodes: KM.ModelPart = model_part.GetSubModelPart("interface_nodes")

        self.assertEqual(smp_interface.NumberOfNodes(), smp_interface_nodes.NumberOfNodes())
        self.assertEqual(smp_interface.NumberOfNodes(), 27)

        self.assertEqual(smp_interface.NumberOfGeometries(), 36)
        self.assertEqual(smp_interface_nodes.NumberOfGeometries(), 0)

        smp_interface_node_ids_unique: Set[int] = set(get_node_ids(smp_interface))
        smp_interface_nodes_node_ids_unique: Set[int] = set(get_node_ids(smp_interface_nodes))

        # make sure that the smp has the nodes of its geometries
        self.assertEqual(smp_interface_node_ids_unique, set(get_geom_node_ids(smp_interface)))

        # check if smps have same nodes
        self.assertEqual(smp_interface_node_ids_unique, smp_interface_nodes_node_ids_unique)

        # check node coordinates
        for node in model_part.Nodes:
            self.assertTrue(0.0 <= node.X <= 200.0)
            self.assertTrue(0.0 <= node.X0 <= 200.0)

            self.assertTrue(0.0 <= node.Y <= 200.0)
            self.assertTrue(0.0 <= node.Y0 <= 200.0)

            self.assertTrue(0.0 <= node.Z <= 200.0)
            self.assertTrue(0.0 <= node.Z0 <= 200.0)

        for node in chain(smp_interface.Nodes, smp_interface_nodes.Nodes):
            self.assertAlmostEqual(node.X, 200.0)
            self.assertAlmostEqual(node.X0, 200.0)

            self.assertTrue(0.0 <= node.Y <= 200.0)
            self.assertTrue(0.0 <= node.Y0 <= 200.0)

            self.assertTrue(0.0 <= node.Z <= 200.0)
            self.assertTrue(0.0 <= node.Z0 <= 200.0)

        # check that the correct geoms (3D triangles) are in the smp
        for geom in smp_interface.Geometries:
            self.assertIsInstance(geom, KM.Triangle3D3)

        # check how many geoms of each type
        exp_geoms: Dict[Any, int] = {
            KM.Tetrahedra3D4: 386,
            KM.Triangle3D3: 210,
            KM.Line3D2: 48
        }
        self.assertEqual(sum(exp_geoms.values()), 644)
        self.assertDictEqual(exp_geoms, get_num_geometries_by_type(model_part))

    def test_cube_with_adjacent_groups(self):
        model: KM.Model = KM.Model()
        model_part: KM.ModelPart = model.CreateModelPart("test")
        KratosMed.MedModelPartIO(GetMedPath("cube_with_adjacent_groups", "cube.med")).ReadModelPart(model_part)

        self._basic_checks(model_part)

        self.assertEqual(model_part.NumberOfNodes(), 50)
        self.assertEqual(model_part.NumberOfGeometries(), 226)

        self.assertEqual(model_part.NumberOfSubModelParts(), 3)
        self.assertTrue(model_part.HasSubModelPart("Face_1"))
        self.assertTrue(model_part.HasSubModelPart("Face_2"))
        self.assertTrue(model_part.HasSubModelPart("Edge_1"))

        smp_face_1: KM.ModelPart = model_part.GetSubModelPart("Face_1")
        smp_face_2: KM.ModelPart = model_part.GetSubModelPart("Face_2")
        smp_edge_1: KM.ModelPart = model_part.GetSubModelPart("Edge_1")

        self.assertEqual(smp_face_1.NumberOfNodes(), 15)
        self.assertEqual(smp_face_2.NumberOfNodes(), 15)
        self.assertEqual(smp_edge_1.NumberOfNodes(), 4)

        self.assertEqual(smp_face_1.NumberOfGeometries(), 16)
        self.assertEqual(smp_face_2.NumberOfGeometries(), 16)
        self.assertEqual(smp_edge_1.NumberOfGeometries(), 3)

        # check node coordinates
        for node in model_part.Nodes:
            self.assertTrue(0.0 <= node.X <= 200.0)
            self.assertTrue(0.0 <= node.X0 <= 200.0)

            self.assertTrue(0.0 <= node.Y <= 200.0)
            self.assertTrue(0.0 <= node.Y0 <= 200.0)

            self.assertTrue(0.0 <= node.Z <= 200.0)
            self.assertTrue(0.0 <= node.Z0 <= 200.0)

        # for node in chain(smp_interface.Nodes, smp_interface_nodes.Nodes):
        #     self.assertAlmostEqual(node.X, 200.0)
        #     self.assertAlmostEqual(node.X0, 200.0)

        #     self.assertTrue(0.0 <= node.Y <= 200.0)
        #     self.assertTrue(0.0 <= node.Y0 <= 200.0)

        #     self.assertTrue(0.0 <= node.Z <= 200.0)
        #     self.assertTrue(0.0 <= node.Z0 <= 200.0)

        # check that the correct geoms (3D triangles) are in the smp
        for geom in chain(smp_face_1.Geometries, smp_face_2.Geometries):
            self.assertIsInstance(geom, KM.Triangle3D3)

        # check that the correct geoms (3D lines) are in the smp
        for geom in smp_edge_1.Geometries:
            self.assertIsInstance(geom, KM.Line3D2)

        # check how many geoms of each type
        exp_geoms: Dict[Any, int] = {
            KM.Tetrahedra3D4: 94,
            KM.Triangle3D3: 96,
            KM.Line3D2: 36
        }
        self.assertEqual(sum(exp_geoms.values()), 226)
        self.assertDictEqual(exp_geoms, get_num_geometries_by_type(model_part))

    def _basic_checks(self, model_part):
        # check no elements or conditions are created
        self.assertEqual(model_part.NumberOfElements(), 0)
        self.assertEqual(model_part.NumberOfConditions(), 0)

        # check increasing node Ids
        for i, node in enumerate(model_part.Nodes):
            self.assertEqual(node.Id, i+1)

        # check geometries have correct Ids
        # Note: Geometries are not ordered, thus cannot check like nodes
        self.assertEqual(set(get_geometry_ids(model_part)), set(range(1, model_part.NumberOfGeometries()+1)))

        # check that the entities are unique in the ModelParts
        self._check_unique_nodes(model_part)
        self._check_unique_geometries(model_part)

        # check each ModelPart has (at least) the nodes of its geometries
        self._check_nodes_geometries(model_part)

    def _check_unique_nodes(self, model_part):
        node_ids: List[int] = get_node_ids(model_part)
        node_ids_unique: Set[int] = set(node_ids)
        self.assertEqual(len(node_ids), len(node_ids_unique))

        for smp in model_part.SubModelParts:
            self._check_unique_nodes(smp)

    def _check_unique_geometries(self, model_part):
        geom_ids: List[int] = get_geometry_ids(model_part)
        geom_ids_unique: Set[int] = set(geom_ids)
        self.assertEqual(len(geom_ids), len(geom_ids_unique))

        for smp in model_part.SubModelParts:
            self._check_unique_geometries(smp)

    def _check_nodes_geometries(self, model_part):
        geom_node_ids: List[int] = get_geom_node_ids(model_part)
        geom_node_ids_unique: Set[int] = set(geom_node_ids)

        node_ids: List[int] = get_node_ids(model_part)

        self.assertTrue(geom_node_ids_unique.issubset(node_ids))

        for smp in model_part.SubModelParts:
            self._check_nodes_geometries(smp)


def get_num_geometries_by_type(model_part: KM.ModelPart) -> Dict[Any, int]:
    geoms_by_type: Dict[Any, int] = {}
    for geom in model_part.Geometries:
        type_geom = type(geom)
        if type_geom not in geoms_by_type:
            geoms_by_type[type_geom] = 0
        geoms_by_type[type_geom] += 1
    return geoms_by_type


def get_node_ids(model_part: KM.ModelPart) -> List[int]:
    return [node.Id for node in model_part.Nodes]


def get_geometry_ids(model_part: KM.ModelPart) -> List[int]:
    return [geom.Id for geom in model_part.Geometries]


def get_geom_node_ids(model_part: KM.ModelPart) -> List[int]:
    node_ids: List[int] = []
    for geom in model_part.Geometries:
        for node in geom:
            node_ids.append(node.Id)
    return node_ids

if __name__ == '__main__':
    KratosUnittest.main()

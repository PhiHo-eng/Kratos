//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <sstream>
#include <vector>
#include <tuple>

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/xml_utilities/xml_element.h"
#include "utilities/xml_utilities/xml_ostream_writer.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(XmlElementGetTagName, KratosCoreFastSuite)
{
    XmlElement element("TestElement");
    KRATOS_CHECK_EQUAL(element.GetTagName(), "TestElement");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementAddAndGetAttributes, KratosCoreFastSuite)
{
    XmlElement element("TestElement");
    element.AddAttribute("Attribute1", "Value1");
    element.AddAttribute("Attribute2", "Value2");

    auto attributes = element.GetAttributes();

    KRATOS_CHECK_EQUAL(attributes[0].first, "Attribute1");
    KRATOS_CHECK_EQUAL(attributes[0].second, "Value1");
    KRATOS_CHECK_EQUAL(attributes[1].first, "Attribute2");
    KRATOS_CHECK_EQUAL(attributes[1].second, "Value2");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementClearAttributes, KratosCoreFastSuite)
{
    XmlElement element("TestElement");
    element.AddAttribute("Attribute1", "Value1");
    element.ClearAttributes();

    auto attributes = element.GetAttributes();
    KRATOS_CHECK_EQUAL(attributes.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementAddElement, KratosCoreFastSuite)
{
    XmlElement element("TestElement");
    auto childElement = Kratos::make_shared<XmlElement>("ChildElement");
    element.AddElement(childElement);

    auto children = element.GetElements();
    KRATOS_CHECK_EQUAL(children[0]->GetTagName(), "ChildElement");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementGetElements, KratosCoreFastSuite)
{
    XmlElement element("TestElement");
    auto childElement1 = Kratos::make_shared<XmlElement>("ChildElement");
    auto childElement2 = Kratos::make_shared<XmlElement>("ChildElement2");
    element.AddElement(childElement1);
    element.AddElement(childElement2);

    auto children = element.GetElements("ChildElement");
    KRATOS_CHECK_EQUAL(children.size(), 1);
    KRATOS_CHECK_EQUAL(children[0]->GetTagName(), "ChildElement");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementClearElements, KratosCoreFastSuite)
{
    XmlElement element("TestElement");
    auto childElement = Kratos::make_shared<XmlElement>("ChildElement");
    element.AddElement(childElement);
    element.ClearElements();

    auto children = element.GetElements();
    KRATOS_CHECK_EQUAL(children.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWrite, KratosCoreFastSuite)
{
    XmlElement element("TestElement");
    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::ASCII, 4);
    element.AddAttribute("Attribute1", "Value1");
    auto childElement = Kratos::make_shared<XmlElement>("ChildElement");
    element.AddElement(childElement);
    element.Write(writer, 1);

    KRATOS_CHECK_EQUAL(ss.str(),
                       "   <TestElement Attribute1=\"Value1\">\n"
                       "      <ChildElement/>\n"
                       "   </TestElement>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteElement, KratosCoreFastSuite)
{
    std::vector<std::pair<std::string, std::string>> attributes = {{"attribute1", "value1"}, {"attribute2", "value2"}};
    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::ASCII, 4);
    writer.WriteElement("TestElement", attributes, 1, true);
    KRATOS_CHECK_EQUAL(ss.str(), "   <TestElement attribute1=\"value1\" attribute2=\"value2\"/>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterCloseElement, KratosCoreFastSuite)
{
    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::ASCII, 4);
    writer.CloseElement("TestElement", 1);
    KRATOS_CHECK_EQUAL(ss.str(), "   </TestElement>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementAsciiChar, KratosCoreFastSuite)
{
    auto char_expression_1 = LiteralFlatExpression<char>::Create(2, {3});
    std::size_t local_index = 0;
    for (auto it = char_expression_1->begin(); it != char_expression_1->end(); ++it) {
        *it = local_index++;
    }

    auto char_expression_2 = LiteralFlatExpression<char>::Create(3, {3});
    local_index = 0;
    for (auto it = char_expression_2->begin(); it != char_expression_2->end(); ++it) {
        *it = local_index++;
    }

    std::vector<Expression::Pointer> expressions = {char_expression_1, char_expression_2};
    XmlElement element("data_1", expressions);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::ASCII, 4);
    element.Write(writer, 1);

    KRATOS_CHECK_EQUAL(ss.str(),
            "   <DataArray type=\"UInt8\" Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\">\n"
            "     0  1  2  3  4  5  0  1  2  3  4  5  6  7  8\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementAsciiInt, KratosCoreFastSuite)
{
    auto char_expression_1 = LiteralFlatExpression<int>::Create(2, {3});
    std::size_t local_index = 0;
    for (auto it = char_expression_1->begin(); it != char_expression_1->end(); ++it) {
        *it = local_index++;
    }

    auto char_expression_2 = LiteralFlatExpression<int>::Create(3, {3});
    local_index = 0;
    for (auto it = char_expression_2->begin(); it != char_expression_2->end(); ++it) {
        *it = local_index++;
    }

    std::vector<Expression::Pointer> expressions = {char_expression_1, char_expression_2};
    XmlElement element("data_1", expressions);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::ASCII, 4);
    element.Write(writer, 1);

    KRATOS_CHECK_EQUAL(ss.str(),
            "   <DataArray type=\"Int32\" Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\">\n"
            "     0  1  2  3  4  5  0  1  2  3  4  5  6  7  8\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementAsciiDouble, KratosCoreFastSuite)
{
    auto char_expression_1 = LiteralFlatExpression<double>::Create(2, {3});
    std::size_t local_index = 0;
    for (auto it = char_expression_1->begin(); it != char_expression_1->end(); ++it) {
        *it = local_index++;
    }

    auto char_expression_2 = LiteralFlatExpression<double>::Create(3, {3});
    local_index = 0;
    for (auto it = char_expression_2->begin(); it != char_expression_2->end(); ++it) {
        *it = local_index++;
    }

    std::vector<Expression::Pointer> expressions = {char_expression_1, char_expression_2};
    XmlElement element("data_1", expressions);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::ASCII, 4);
    element.Write(writer, 1);

    KRATOS_CHECK_EQUAL(ss.str(),
            "   <DataArray type=\"Float64\" Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\">\n"
            "     0.0000e+00  1.0000e+00  2.0000e+00  3.0000e+00  4.0000e+00  5.0000e+00  0.0000e+00  1.0000e+00  2.0000e+00  3.0000e+00  4.0000e+00  5.0000e+00  6.0000e+00  7.0000e+00  8.0000e+00\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementAsciiMixed, KratosCoreFastSuite)
{
    auto char_expression_1 = LiteralFlatExpression<char>::Create(2, {3});
    std::size_t local_index = 0;
    for (auto it = char_expression_1->begin(); it != char_expression_1->end(); ++it) {
        *it = local_index++;
    }

    auto char_expression_2 = LiteralFlatExpression<double>::Create(3, {3});
    local_index = 0;
    for (auto it = char_expression_2->begin(); it != char_expression_2->end(); ++it) {
        *it = local_index++;
    }

    std::vector<Expression::Pointer> expressions = {char_expression_1, char_expression_2};
    XmlElement element("data_1", expressions);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::ASCII, 4);
    element.Write(writer, 1);

    KRATOS_CHECK_EQUAL(ss.str(),
            "   <DataArray type=\"Float64\" Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\">\n"
            "     0.0000e+00  1.0000e+00  2.0000e+00  3.0000e+00  4.0000e+00  5.0000e+00  0.0000e+00  1.0000e+00  2.0000e+00  3.0000e+00  4.0000e+00  5.0000e+00  6.0000e+00  7.0000e+00  8.0000e+00\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementBinaryChar, KratosCoreFastSuite)
{
    auto char_expression_1 = LiteralFlatExpression<char>::Create(2, {3});
    std::size_t local_index = 0;
    for (auto it = char_expression_1->begin(); it != char_expression_1->end(); ++it) {
        *it = local_index++;
    }

    auto char_expression_2 = LiteralFlatExpression<char>::Create(3, {3});
    local_index = 0;
    for (auto it = char_expression_2->begin(); it != char_expression_2->end(); ++it) {
        *it = local_index++;
    }

    std::vector<Expression::Pointer> expressions = {char_expression_1, char_expression_2};
    XmlElement element("data_1", expressions);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::BINARY, 4);
    element.Write(writer, 1);

        KRATOS_CHECK_EQUAL(ss.str(),
            "   <DataArray type=\"UInt8\" Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\">\n"
            "     DwAAAAABAgMEBQABAgMEBQYHCA==\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementBinaryInt, KratosCoreFastSuite)
{
    auto char_expression_1 = LiteralFlatExpression<int>::Create(2, {3});
    std::size_t local_index = 0;
    for (auto it = char_expression_1->begin(); it != char_expression_1->end(); ++it) {
        *it = local_index++;
    }

    auto char_expression_2 = LiteralFlatExpression<int>::Create(3, {3});
    local_index = 0;
    for (auto it = char_expression_2->begin(); it != char_expression_2->end(); ++it) {
        *it = local_index++;
    }

    std::vector<Expression::Pointer> expressions = {char_expression_1, char_expression_2};
    XmlElement element("data_1", expressions);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::BINARY, 4);
    element.Write(writer, 1);

    KRATOS_CHECK_EQUAL(ss.str(),
            "   <DataArray type=\"Int32\" Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\">\n"
            "     PAAAAAAAAAABAAAAAgAAAAMAAAAEAAAABQAAAAAAAAABAAAAAgAAAAMAAAAEAAAABQAAAAYAAAAHAAAACAAAAA==\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementBinaryDouble, KratosCoreFastSuite)
{
    auto char_expression_1 = LiteralFlatExpression<double>::Create(2, {3});
    std::size_t local_index = 0;
    for (auto it = char_expression_1->begin(); it != char_expression_1->end(); ++it) {
        *it = local_index++;
    }

    auto char_expression_2 = LiteralFlatExpression<double>::Create(3, {3});
    local_index = 0;
    for (auto it = char_expression_2->begin(); it != char_expression_2->end(); ++it) {
        *it = local_index++;
    }

    std::vector<Expression::Pointer> expressions = {char_expression_1, char_expression_2};
    XmlElement element("data_1", expressions);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::BINARY, 4);
    element.Write(writer, 1);

    KRATOS_CHECK_EQUAL(ss.str(),
            "   <DataArray type=\"Float64\" Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\">\n"
            "     eAAAAAAAAAAAAAAAAAAAAAAA8D8AAAAAAAAAQAAAAAAAAAhAAAAAAAAAEEAAAAAAAAAUQAAAAAAAAAAAAAAAAAAA8D8AAAAAAAAAQAAAAAAAAAhAAAAAAAAAEEAAAAAAAAAUQAAAAAAAABhAAAAAAAAAHEAAAAAAAAAgQA==\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementBinaryMixed, KratosCoreFastSuite)
{
    auto char_expression_1 = LiteralFlatExpression<char>::Create(2, {3});
    std::size_t local_index = 0;
    for (auto it = char_expression_1->begin(); it != char_expression_1->end(); ++it) {
        *it = local_index++;
    }

    auto char_expression_2 = LiteralFlatExpression<double>::Create(3, {3});
    local_index = 0;
    for (auto it = char_expression_2->begin(); it != char_expression_2->end(); ++it) {
        *it = local_index++;
    }

    std::vector<Expression::Pointer> expressions = {char_expression_1, char_expression_2};
    XmlElement element("data_1", expressions);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    XmlOStreamWriter writer(ss, XmlOStreamWriter::WriterFormat::BINARY, 4);
    element.Write(writer, 1);

    KRATOS_CHECK_EQUAL(ss.str(),
            "   <DataArray type=\"Float64\" Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\">\n"
            "     eAAAAAAAAAAAAAAAAAAAAAAA8D8AAAAAAAAAQAAAAAAAAAhAAAAAAAAAEEAAAAAAAAAUQAAAAAAAAAAAAAAAAAAA8D8AAAAAAAAAQAAAAAAAAAhAAAAAAAAAEEAAAAAAAAAUQAAAAAAAABhAAAAAAAAAHEAAAAAAAAAgQA==\n"
            "   </DataArray>\n");
}

} // namespace Kratos::Testing
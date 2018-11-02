//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#if !defined(KRATOS_FILE_SERIALIZER_H_INCLUDED )
#define  KRATOS_FILE_SERIALIZER_H_INCLUDED

// System includes
#include <string>
#include <cstring>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{
    //This is a Helper object to provide a simpler interface for serialization to a file
    class KRATOS_API(KRATOS_CORE) FileSerializer : public Serializer 
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(FileSerializer); 

            ///this constructor simply wraps the standard Serializer and defines output to basic_iostream
            ///@param rTrace type of serialization to be employed
            FileSerializer(std::string const& Filename, Serializer::TraceType const& rTrace=SERIALIZER_NO_TRACE) 
                : Serializer(nullptr, rTrace)
            {
                std::fstream* p_file = new std::fstream(std::string(Filename+".rest").c_str(), std::ios::binary|std::ios::in|std::ios::out);
                if(!(*p_file))
                {
                    delete p_file;
                    p_file = new std::fstream(std::string(Filename+".rest").c_str(), std::ios::binary|std::ios::out);
                }
                SetBuffer( p_file );
                KRATOS_ERROR_IF_NOT(*pGetBuffer()) << "Error opening input file : "
                                            << std::string(Filename+".rest") << std::endl;
            }

            virtual ~FileSerializer(){}

        private:

            /// Assignment operator.
            FileSerializer& operator=(FileSerializer const& rOther) = delete;

            /// Copy constructor.
            FileSerializer(FileSerializer const& rOther) = delete;
    };
}

#endif // KRATOS_FILE_SERIALIZER_H_INCLUDED  defined

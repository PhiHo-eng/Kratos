//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen, https://github.com/matekelemen
//

#ifndef KRATOS_PATTERN_UTILITY_H_INCLUDED
#define KRATOS_PATTERN_UTILITY_H_INCLUDED

// External includes
#include "ghc/filesystem.hpp" // TODO change to std::filesystem when migrating to C++17

// Core includes
#include "includes/kratos_export_api.h"
#include "includes/shared_pointers.h"
#include "includes/model_part.h"

// STL includes
#include <regex>
#include <vector>
#include <map>


namespace Kratos
{


///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/** A class for interfacing placeholders and regular expressions.
 *
 *  @note placeholders should be separated by literals, otherwise
 *  regex will probably not capture them as you'd expect.
 *  BAD example:
 *      pattern:    "<placeholder_1><placeholder_2>"
 *      string:     "abcdefg"
 *      result:
 *          "<placeholder_1>" : "abcdef"
 *          "<placeholder_2>" : "g"
 *  CORRECT example:
 *      pattern:    "<placeholder_1>d<placeholder_2>"
 *      string:     "abcdefg"
 *      result:
 *          "<placeholder_1>" : "abc"
 *          "<placeholder_2>" : "efg"
 *
 *  @note placeholders not present in the pattern are discarded,
 *  and won't be keys in @ref{PlaceholderPattern::Match}.
 */
class KRATOS_API(HDF5Application) PlaceholderPattern
{
public:
    ///@name Type Definitions
    ///@{

    using PlaceholderMap = std::map<std::string,std::string>;

    using PlaceholderGroupMap = std::map<std::string,std::vector<std::size_t>>;

    using MatchType = std::map<std::string,std::vector<std::string>>;

    using PathType = ghc::filesystem::path;

    KRATOS_CLASS_POINTER_DEFINITION(PlaceholderPattern);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    PlaceholderPattern() = default;

    /** Constructor
     *  @param rPattern pattern string with placeholders
     *  @param rPlaceholders pairs of placeholders and their corresponding regex strings
     *                       Example: {{"<name>", ".+"}, {"<identifier>", "[0-9]+"}}
     *
     *  @note the corresponding regexes must be bare, not containing groups or position
     *  constraints (such as line begin or end modifiers).
     */
    PlaceholderPattern(const std::string& rPattern,
                       const PlaceholderMap& rPlaceholderMap);

    /// Move constructor
    PlaceholderPattern(PlaceholderPattern&& rOther) = default;

    /// Copy constructor
    PlaceholderPattern(const PlaceholderPattern& rOther) = default;

    virtual ~PlaceholderPattern() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Move assignment operator
    PlaceholderPattern& operator=(PlaceholderPattern&& rOther) = default;

    /// Copy assignment operator
    PlaceholderPattern& operator=(const PlaceholderPattern& rOther) = default;

    ///@}
    ///@name Operations
    ///@{

    /// Check whether a string satisfies the pattern
    bool IsAMatch(const std::string& rString) const;

    /** Find all placeholders' values in the input string.
     *  @param rString string that matches the input pattern
     *  @return a map associating a vector of strings, i.e. the values
     *  of placeholders in the input string, to the placeholders
     *
     *  @note the returned placeholder values appear in the same order
     *  they do in the input pattern.
     */
    MatchType Match(const std::string& rString) const;

    /**
     *  Return a copy of the pattern that has its placeholders replaced
     *  with the corresponding values specified in the input map.
     *
     *  @param rPlaceholderValueMap string - string map associating values to placeholders
     *                              {"palceholder" : "placeholder_value"}
     */
    std::string Apply(const PlaceholderMap& rPlaceholderValueMap) const;

    /** Collect all file/directory paths that match the pattern.
     *
     *  @note the search begins from the filesystem root if the pattern is an absolute path,
     *  otherwise it begins from the cwd.
     */
    std::vector<PathType> Glob() const;

    ///@}
    ///@name Inquiry
    ///@{

    /// Get the regex for the input pattern
    const std::regex& GetRegex() const;

    /// Get the string representation of the regex
    const std::string& GetRegexString() const;

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mPattern;

    // Placeholders and their assigned group indices in the pattern
    PlaceholderGroupMap mPlaceholderGroupMap;

    std::string mRegexString;

    std::regex mRegex;

    ///@}
    ///@name Private Operations
    ///@{

    /// Escape characters in the input that interfere with regex syntax
    static std::string FormatRegexLiteral(const std::string& rLiteral);

    ///@}
}; // class PlaceholderPattern

///@} // Kratos Classes
///@} addtogroup
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/** A class for working with formatted strings related to @ref{ModelPart}s.
 *
 *  Operations on strings with the following placeholders are supported:
 *      <model_part_name>
 *      <step>
 *      <time>
 *  See @ref{PlaceholderPattern} for supported functionalities. Other
 *  placeholders can be added at compile time by tweaking the construction
 *  of the static member @ref{ModelPartPattern::mModelpartPlaceholderMap}.
 */
class KRATOS_API(HDF5Application) ModelPartPattern : public PlaceholderPattern
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ModelPartPattern);

    ///@}
    ///@name Life Cycle
    ///@{

    ModelPartPattern() = default;

    ModelPartPattern(const std::string& rPattern);

    ModelPartPattern(ModelPartPattern&& rOther) = default;

    ModelPartPattern(const ModelPartPattern& rOther) = default;

    ///@}
    ///@name Operators
    ///@{

    using PlaceholderPattern::operator=;

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static const PlaceholderPattern::PlaceholderMap mModelPartPlaceholderMap;

    ///@}
}; // class ModelPartPattern

///@} // Kratos Classes
///@} addtogroup


} // namespace Kratos


#endif

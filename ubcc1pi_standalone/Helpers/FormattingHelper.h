/**
 *  @file  ubcc1pi_standalone/Helpers/FormattingHelper.h
 *
 *  @brief The header file for the formatting helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_FORMATTING_HELPER
#define UBCC1PI_STANDALONE_HELPERS_FORMATTING_HELPER

#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>

namespace ubcc1pi
{

/**
 *  @brief  The formatting helper class
 */
class FormattingHelper
{
    public:
        /**
         *  @brief  The table class
         */
        class Table
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  headers the headers of the table columns, use empty string for spacer columns
                 */
                Table(const std::vector<std::string> &headers);

                /**
                 *  @brief  Get the number of rows in the table
                 *
                 *  @return the number of rows
                 */
                unsigned int GetNumberOfRows() const;

                /**
                 *  @brief  Print the table in markdown syntax
                 */
                void Print() const;

                /**
                 *  @brief  Set the entry for the given header and row number
                 *
                 *  @tparam T the type of entry to add
                 *  @param  header the header to set
                 *  @param  row the row to set
                 *  @param  value the value to set
                 */
                template <typename T>
                    void SetEntry(const std::string &header, const unsigned int row, const T &value);

                /**
                 *  @brief  Set the entry for the given header on the last row
                 *
                 *  @tparam T the type of entry to add
                 *  @param  header the header to set
                 *  @param  value the value to set
                 */
                template <typename T>
                    void SetEntry(const std::string &header, const T &value);

                /**
                 *  @brief  Add a new empty row
                 */
                void AddEmptyRow();

            private:

                /**
                 *  @brief  Get the widths of a column for printing
                 *
                 *  @param  the header of the column
                 *
                 *  @return the column width
                 */
                unsigned int GetColumnWidth(const std::string &header) const;

                std::vector<std::string>                                     m_headers;      ///< The column headers
                std::vector<std::unordered_map< std::string, std::string > > m_entries;      ///< The table entries [row][column]
        };

        /**
         *  @brief  Print a line of '-' characters
         */
        static void PrintLine();
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

FormattingHelper::Table::Table(const std::vector<std::string> &headers) :
    m_headers(headers)
{
    for (const auto &header : m_headers)
    {
        if (header.empty())
            continue;

        const auto nEntries = std::count(m_headers.begin(), m_headers.end(), header);
        if (nEntries != 1)
            throw std::invalid_argument("FormattingHelper::Table::Table - Repeated header: \"" + header + "\"");
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
            
unsigned int FormattingHelper::Table::GetNumberOfRows() const
{
    return m_entries.size();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int FormattingHelper::Table::GetColumnWidth(const std::string &header) const
{
    if (header.empty())
        return 0u;
    
    auto width = header.size();
    for (const auto &entriesInRow : m_entries)
    {
        const auto &entry = entriesInRow.at(header);
        width = std::max(width, entry.size());
    }
    
    return width;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FormattingHelper::Table::Print() const
{
    std::vector<unsigned int> widths;

    // Print the header row
    for (const auto &header : m_headers)
    {
        const auto width = this->GetColumnWidth(header);
        widths.push_back(width);

        std::cout << "|";

        if (header.empty())
            continue;

        std::cout << " " << std::setw(width) << std::left << header << " ";
    }
    std::cout << "|" << std::endl;
    
    // Print the horizontal line
    for (unsigned int iHeader = 0; iHeader < m_headers.size(); ++iHeader)
    {
        const auto &header = m_headers.at(iHeader);

        std::cout << "|";

        if (header.empty())
            continue;
        
        const auto width = widths.at(iHeader);

        std::cout << std::string(width + 2, '-');
    }
    std::cout << "|" << std::endl;

    // Print the entries
    for (const auto &entriesInRow : m_entries)
    {
        for (unsigned int iHeader = 0; iHeader < m_headers.size(); ++iHeader)
        {
            const auto &header = m_headers.at(iHeader);

            std::cout << "|";

            if (header.empty())
                continue;
            
            const auto width = widths.at(iHeader);
            const auto &entry = entriesInRow.at(header);

            std::cout << " " << std::setw(width) << std::left << entry << " ";
        }
        std::cout << "|" << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void FormattingHelper::Table::SetEntry(const std::string &header, const T &value)
{
    if (m_entries.empty())
        throw std::invalid_argument("FormattingHelper::Table::SetEntry - Table has now rows");

    this->SetEntry(header, this->GetNumberOfRows() - 1, value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void FormattingHelper::Table::SetEntry(const std::string &header, const unsigned int row, const T &value)
{
    if (header.empty())
        throw std::invalid_argument("FormattingHelper::Table::SetEntry - Input header is empty string");

    if (row >= m_entries.size()) 
        throw std::invalid_argument("FormattingHelper::Table::SetEntry - Input row is out of bounds");

    // Find this header
    auto &entriesInRow = m_entries.at(row);
    auto entryIter = entriesInRow.find(header);
    if (entryIter == entriesInRow.end())
        throw std::invalid_argument("FormattingHelper::Table::SetEntry - Invalid input header: \"" + header + "\"");

    auto &entry = entryIter->second;

    // Set the entry to the value
    stringstream ss;
    ss.precision(5);
    ss << value;
    entry = ss.str();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
            
void FormattingHelper::Table::AddEmptyRow()
{
    std::unordered_map<std::string, std::string> row;
    for (const auto &header : m_headers)
    {
        if (header.empty())
            continue;

        row.emplace(header, "");
    }

    m_entries.push_back(row);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void FormattingHelper::PrintLine()
{
    std::cout << std::string(87, '-') << std::endl;
}

} // namespace ubcc1pi

#endif

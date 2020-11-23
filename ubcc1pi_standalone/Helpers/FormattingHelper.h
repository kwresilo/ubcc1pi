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
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#include <TH1F.h>
#include <TH2F.h>

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
                 *  @brief  Write the table in markdown syntax to an output file
                 *
                 *  @param  fileName the output file name
                 *  @param  alsoPrint if we should also print the table to the terminal
                 */
                void WriteToFile(const std::string &fileName, const bool alsoPrint = true) const;

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
                 *  @param  header the header of the column
                 *
                 *  @return the column width
                 */
                unsigned int GetColumnWidth(const std::string &header) const;

                std::vector<std::string>                                     m_headers;      ///< The column headers
                std::vector<std::unordered_map< std::string, std::string > > m_entries;      ///< The table entries [row][column]
        };

        /**
         *  @brief  Get the string of the form "3.14 +- 0.03" showing a value with the error
         *
         *  @param  value the value
         *  @param  uncertainty the uncertainty
         *
         *  @return the value with error
         */
        static std::string GetValueWithError(const float &value, const float &uncertainty);

        /**
         *  @brief  Save a 1D histogram as a table
         *
         *  @param  pHist the input histogram
         *  @param  fileName the output file name
         *  @param  alsoPrint if we should also print the table to the terminal
         */
        static void SaveHistAsTable(const std::shared_ptr<TH1F> &pHist, const std::string &fileName, const bool alsoPrint = true);

        /**
         *  @brief  Save a 2D histogram as a table
         *
         *  @param  pHist the input histogram
         *  @param  fileName the output file name
         *  @param  alsoPrint if we should also print the table to the terminal
         */
        static void SaveHistAsTable(const std::shared_ptr<TH2F> &pHist, const std::string &fileName, const bool alsoPrint = true);

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

void FormattingHelper::Table::WriteToFile(const std::string &fileName, const bool alsoPrint) const
{
    fstream file;
    file.open(fileName, ios::out);

    // Get the stream buffers for cout and the file
    std::streambuf* stream_buffer_cout = std::cout.rdbuf();
    std::streambuf* stream_buffer_file = file.rdbuf();

    // Redirect cout to the file
    cout.rdbuf(stream_buffer_file);
    this->Print();

    // Redirect cout back to the terminal
    cout.rdbuf(stream_buffer_cout);
    file.close();

    // Also print to the terminal if requested
    if (alsoPrint)
        this->Print();
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
// -----------------------------------------------------------------------------------------------------------------------------------------

std::string FormattingHelper::GetValueWithError(const float &value, const float &uncertainty)
{
    stringstream ss;

    if (uncertainty < std::numeric_limits<float>::epsilon())
    {
        // Just convert the floats to strings
        ss << value << " +- " << uncertainty;
    }
    else
    {
        // Determine the precision we should use
        const auto nDecimalPlaces = std::floor(std::log10(uncertainty)) - 1;
        const auto powerOfTen = std::pow(10.f, nDecimalPlaces);

        // Round to that precision
        //const auto valueRounded = powerOfTen * std::round(value / powerOfTen);
        const auto uncertaintyRounded = powerOfTen * std::round(uncertainty / powerOfTen);

        // Convert the floats to strings
        ss << value << " +- " << uncertaintyRounded;
    }

    return ss.str();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FormattingHelper::PrintLine()
{
    std::cout << std::string(87, '-') << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FormattingHelper::SaveHistAsTable(const std::shared_ptr<TH1F> &pHist, const std::string &fileName, const bool alsoPrint)
{
    FormattingHelper::Table table({"Bin", "Lower", "Upper", "Width", "", "Value"});
    for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pHist->GetNbinsX()); ++iBin)
    {
        table.AddEmptyRow();
        table.SetEntry("Bin", iBin);

        const auto lower = pHist->GetBinLowEdge(iBin);
        table.SetEntry("Lower", lower);

        const auto width = pHist->GetBinWidth(iBin);
        table.SetEntry("Width", width);

        const auto upper = lower + width;
        table.SetEntry("Upper", upper);

        const auto value = pHist->GetBinContent(iBin);
        table.SetEntry("Value", value);
    }

    table.WriteToFile(fileName, alsoPrint);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FormattingHelper::SaveHistAsTable(const std::shared_ptr<TH2F> &pHist, const std::string &fileName, const bool alsoPrint)
{
    FormattingHelper::Table table({"Bin X", "Bin Y", "Lower", "Upper", "Width", "", "Value"});

    for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pHist->GetNbinsX()); ++iBin)
    {
        table.AddEmptyRow();
        table.SetEntry("Bin X", iBin);

        const auto xAxis = pHist->GetXaxis();
        const auto xLower = xAxis->GetBinLowEdge(iBin);
        table.SetEntry("Lower", xLower);

        const auto xWidth = xAxis->GetBinWidth(iBin);
        table.SetEntry("Width", xWidth);

        const auto xUpper = xLower + xWidth;
        table.SetEntry("Upper", xUpper);

        for (unsigned int jBin = 1; jBin <= static_cast<unsigned int>(pHist->GetNbinsY()); ++jBin)
        {
            table.AddEmptyRow();
            table.SetEntry("Bin Y", jBin);

            const auto yAxis = pHist->GetYaxis();
            const auto yLower = yAxis->GetBinLowEdge(jBin);
            table.SetEntry("Lower", yLower);

            const auto yWidth = yAxis->GetBinWidth(jBin);
            table.SetEntry("Width", yWidth);

            const auto yUpper = yLower + yWidth;
            table.SetEntry("Upper", yUpper);

            const auto value = pHist->GetBinContent(iBin, jBin);
            table.SetEntry("Value", value);
        }
    }

    table.WriteToFile(fileName, alsoPrint);
}

} // namespace ubcc1pi

#endif

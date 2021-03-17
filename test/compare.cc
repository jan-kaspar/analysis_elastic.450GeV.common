#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"

#include <string>
#include <tuple>
#include <memory>

using namespace std;

//----------------------------------------------------------------------------------------------------

std::tuple<double, bool, bool> CompareValues(const double v1, const double v2)
{
    const double rel_diff = (v1 == 0. && v2 == 0.) ? 0. : (v2 - v1) / (v2 + v1) * 2.;
    const bool diff_minor = (fabs(rel_diff) > 1E-6);
    const bool diff_major = (fabs(rel_diff) > 1E-3);

    return {rel_diff, diff_major, diff_minor};
}

//----------------------------------------------------------------------------------------------------

unique_ptr<TCanvas> NewDefaultCanvas()
{
    return make_unique<TCanvas>("cmp", "", 2048, 1536);
}

//----------------------------------------------------------------------------------------------------

const string title = "(test) blue solid, (ref) red dashed";

void SetStyle(TObject *o, unsigned int style)
{
    int color = 0;
    int lineStyle = 0;
    int width = 2;

    if (style == 1)
    {
        color = 4;
        lineStyle = 1;
    }

    if (style == 2)
    {
        color = 2;
        lineStyle = 2;
    }

    if (o->InheritsFrom("TGraph"))
    {
        TGraph *g = (TGraph *) o;
        g->SetLineStyle(lineStyle);
        g->SetLineColor(color);
        g->SetMarkerColor(color);
        g->SetLineWidth(width);

        for (auto *o : * g->GetListOfFunctions())
           SetStyle(o, style) ;
    }

    if (o->InheritsFrom("TH1"))
    {
        TH1 *h = (TH1 *) o;
        h->SetLineStyle(lineStyle);
        h->SetLineColor(color);
        h->SetLineWidth(width);

        for (auto *o : * h->GetListOfFunctions())
           SetStyle(o, style) ;
    }

    if (o->InheritsFrom("TF1"))
    {
        TF1 *f = (TF1 *) o;
        f->SetLineStyle(lineStyle);
        f->SetLineColor(color);
    }
}

//----------------------------------------------------------------------------------------------------

int CompareGraphsImp(const TGraph *g1, const TGraph *g2)
{
    if (g1->GetN() != g2->GetN())
    {
        printf("* number of points is different\n");
        return 1;
    }

    bool diff_major = false;
    bool diff_minor = false;

    for (int i = 0; i < g1->GetN(); ++i)
    {
        const double x_1 = g1->GetX()[i];
        const double y_1 = g1->GetY()[i];

        const double x_2 = g2->GetX()[i];
        const double y_2 = g2->GetY()[i];

        auto [x_rel_diff, x_diff_major, x_diff_minor] = CompareValues(x_1, x_2);
        diff_major |= x_diff_major;
        diff_minor |= x_diff_minor;
        if (x_diff_minor != 0)
            printf("* point %i: difference in x: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", i, x_1, x_2, x_rel_diff);

        auto [y_rel_diff, y_diff_major, y_diff_minor] = CompareValues(y_1, y_2);
        diff_major |= y_diff_major;
        diff_minor |= y_diff_minor;
        if (y_diff_minor != 0)
            printf("* point %i: difference in x: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", i, y_1, y_2, y_rel_diff);
    }

    if (diff_major)
        return 2;

    if (diff_minor)
        return 3;

    return 0;
}

//----------------------------------------------------------------------------------------------------

int CompareGraphs(TGraph *g1, TGraph *g2, const string& plot_fn)
{
    const int r = CompareGraphsImp(g1, g2);

    if (r != 0 && !plot_fn.empty())
    {
        unique_ptr<TCanvas> c = NewDefaultCanvas();

        SetStyle(g1, 1);
        g1->SetTitle(title.c_str());
        g1->Draw("alp");

        SetStyle(g2, 2);
        g2->Draw("lp");

        c->SaveAs(plot_fn.c_str());
    }

    return r;
}

//----------------------------------------------------------------------------------------------------

int CompareHistograms1DImp(const TH1 *h1, const TH1 *h2)
{
    if (h1->GetNbinsX() != h2->GetNbinsX())
    {
        printf("* number of bins is different\n");
        return 1;
    }

    bool diff_major = false;
    bool diff_minor = false;

    for (int bi = 1; bi <= h1->GetNbinsX(); ++bi)
    {
        const double c_1 = h1->GetBinCenter(bi);
        const double w_1 = h1->GetBinWidth(bi);
        const double v_1 = h1->GetBinContent(bi);

        const double c_2 = h2->GetBinCenter(bi);
        const double w_2 = h2->GetBinWidth(bi);
        const double v_2 = h2->GetBinContent(bi);

        auto [c_rel_diff, c_diff_major, c_diff_minor] = CompareValues(c_1, c_2);
        diff_major |= c_diff_major;
        diff_minor |= c_diff_minor;
        if (c_diff_minor != 0)
            printf("* bin %i: difference in bin center: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", bi, c_1, c_2, c_rel_diff);

        auto [w_rel_diff, w_diff_major, w_diff_minor] = CompareValues(w_1, w_2);
        diff_major |= w_diff_major;
        diff_minor |= w_diff_minor;
        if (w_diff_minor != 0)
            printf("* bin %i: difference in bin width: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", bi, w_1, w_2, w_rel_diff);

        auto [v_rel_diff, v_diff_major, v_diff_minor] = CompareValues(v_1, v_2);
        diff_major |= v_diff_major;
        diff_minor |= v_diff_minor;
        if (v_diff_minor != 0)
            printf("* bin %i: difference in bin content: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", bi, v_1, v_2, v_rel_diff);
    }

    if (diff_major)
        return 2;

    if (diff_minor)
        return 3;

    return 0;
}

//----------------------------------------------------------------------------------------------------

int CompareHistograms1D(TH1 *h1, TH1 *h2, const string &plot_fn)
{
    const int r = CompareHistograms1DImp(h1, h2);

    if (r != 0 && !plot_fn.empty())
    {
        unique_ptr<TCanvas> c = NewDefaultCanvas();

        SetStyle(h1, 1);
        h1->SetTitle(title.c_str());
        h1->Draw("");

        SetStyle(h2, 2);
        h2->Draw("same");

        c->SaveAs(plot_fn.c_str());
    }

    return r;
}

//----------------------------------------------------------------------------------------------------

int CompareHistograms2DImp(const TH2 *h_1, const TH2 *h_2)
{
    if (h_1->GetNbinsX() != h_2->GetNbinsX() || h_1->GetNbinsY() != h_2->GetNbinsY())
    {
        printf("* number of bins is different\n");
        return 1;
    }

    bool diff_major = false;
    bool diff_minor = false;

    // check x bins
    for (int bi = 1; bi <= h_1->GetNbinsX(); ++bi)
    {
        const double c_1 = h_1->GetXaxis()->GetBinCenter(bi);
        const double w_1 = h_1->GetXaxis()->GetBinWidth(bi);

        const double c_2 = h_2->GetXaxis()->GetBinCenter(bi);
        const double w_2 = h_2->GetXaxis()->GetBinWidth(bi);

        auto [c_rel_diff, c_diff_major, c_diff_minor] = CompareValues(c_1, c_2);
        diff_major |= c_diff_major;
        diff_minor |= c_diff_minor;
        if (c_diff_minor != 0)
            printf("* bin x %i: difference in bin center: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", bi, c_1, c_2, c_rel_diff);

        auto [w_rel_diff, w_diff_major, w_diff_minor] = CompareValues(w_1, w_2);
        diff_major |= w_diff_major;
        diff_minor |= w_diff_minor;
        if (w_diff_minor != 0)
            printf("* bin x %i: difference in bin width: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", bi, w_1, w_2, w_rel_diff);
    }

    // check y bins
    for (int bi = 1; bi <= h_1->GetNbinsY(); ++bi)
    {
        const double c_1 = h_1->GetYaxis()->GetBinCenter(bi);
        const double w_1 = h_1->GetYaxis()->GetBinWidth(bi);

        const double c_2 = h_2->GetYaxis()->GetBinCenter(bi);
        const double w_2 = h_2->GetYaxis()->GetBinWidth(bi);

        auto [c_rel_diff, c_diff_major, c_diff_minor] = CompareValues(c_1, c_2);
        diff_major |= c_diff_major;
        diff_minor |= c_diff_minor;
        if (c_diff_minor != 0)
            printf("* bin y %i: difference in bin center: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", bi, c_1, c_2, c_rel_diff);

        auto [w_rel_diff, w_diff_major, w_diff_minor] = CompareValues(w_1, w_2);
        diff_major |= w_diff_major;
        diff_minor |= w_diff_minor;
        if (w_diff_minor != 0)
            printf("* bin y %i: difference in bin width: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", bi, w_1, w_2, w_rel_diff);
    }

    // check bin content
    for (int bi_x = 1; bi_x <= h_1->GetNbinsX(); ++bi_x)
    {
        for (int bi_y = 1; bi_y <= h_1->GetNbinsY(); ++bi_y)
        {
            const int bi = h_1->GetBin(bi_x, bi_y);
            const double v_1 = h_1->GetBinContent(bi);
            const double v_2 = h_2->GetBinContent(bi);

            auto [v_rel_diff, v_diff_major, v_diff_minor] = CompareValues(v_1, v_2);
            diff_major |= v_diff_major;
            diff_minor |= v_diff_minor;
            if (v_diff_minor != 0)
                printf("* bin x %i, bin y %i: difference in bin content: (test) = %.3E, (ref) = %.3E --> rel. diff. = %.1E\n", bi_x, bi_y, v_1, v_2, v_rel_diff);
        }
    }

    if (diff_major)
        return 2;

    if (diff_minor)
        return 3;

    return 0;
}

//----------------------------------------------------------------------------------------------------

int CompareHistograms2D(TH2 *h_1, TH2 *h_2, const string &plot_fn)
{
    const int r = CompareHistograms2DImp(h_1, h_2);

    if (r != 0 && !plot_fn.empty())
    {
        unique_ptr<TCanvas> c = NewDefaultCanvas();
        c->Divide(2, 1);

        c->cd(1);
        h_1->SetTitle("(test)");
        h_1->Draw("colz");

        c->cd(2);
        h_2->SetTitle("(ref)");
        h_2->Draw("colz");

        c->SaveAs(plot_fn.c_str());
    }

    return r;
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
    printf("USAGE: compare <test file> <test object> <ref file> <ref object> [plot file]\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
    if (argc < 5)
    {
        printf("ERROR: insufficient input.");
        PrintUsage();
        return 1;
    }

    string n_file1 = argv[1];
    string n_obj1 = argv[2];
    string n_file2 = argv[3];
    string n_obj2 = argv[4];

    string plot_file = "";
    if (argc >= 6)
        plot_file = argv[5];

    printf("comparing (test) to (ref):\n");
    printf("    (test): %s from %s\n", n_obj1.c_str(), n_file1.c_str());
    printf("    (ref) : %s from %s\n", n_obj2.c_str(), n_file2.c_str());

    // load data
    TFile *file1 = TFile::Open(n_file1.c_str());
    if (!file1)
    {
        printf("ERROR: cannot open file '%s'.\n", n_file1.c_str());
        return 10;
    }

    TObject *obj1 = file1->Get(n_obj1.c_str());
    if (!obj1)
    {
        printf("ERROR: cannot load object '%s' from file1.\n", n_obj1.c_str());
        return 11;
    }

    TFile *file2 = TFile::Open(n_file2.c_str());
    if (!file2)
    {
        printf("ERROR: cannot open file '%s'.\n", n_file2.c_str());
        return 12;
    }

    TObject *obj2 = file2->Get(n_obj2.c_str());
    if (!obj2)
    {
        printf("ERROR: cannot load object '%s' from file2.\n", n_obj2.c_str());
        return 13;
    }

    // do comparison
    const bool obj1_graph = obj1->InheritsFrom("TGraph");
    const bool obj1_histogram = obj1->InheritsFrom("TH1");
    const bool obj1_histogram_2d = obj1->InheritsFrom("TH2");

    const bool obj2_graph = obj2->InheritsFrom("TGraph");
    const bool obj2_histogram = obj2->InheritsFrom("TH1");
    const bool obj2_histogram_2d = obj2->InheritsFrom("TH2");

    if (obj1_graph && obj2_graph)
        return CompareGraphs((TGraph *) obj1, (TGraph *) obj2, plot_file);

    if (obj1_histogram_2d && obj2_histogram_2d)
        return CompareHistograms2D((TH2 *) obj1, (TH2 *) obj2, plot_file);

    if (obj1_histogram && obj2_histogram)
        return CompareHistograms1D((TH1 *) obj1, (TH1 *) obj2, plot_file);

    printf("ERROR: don't know how to compare objects of types %s (test) and %s (ref).\n", obj1->ClassName(), obj2->ClassName());
    return 10;
}

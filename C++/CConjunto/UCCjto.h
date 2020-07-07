//---------------------------------------------------------------------------

#ifndef UCCjtoH
#define UCCjtoH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Buttons.hpp>
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
        TRadioGroup *Cjto;
        TGroupBox *GroupBox1;
        TLabel *Label1;
        TEdit *Edit1;
        TButton *Button1;
        TGroupBox *GroupBox3;
        TLabel *Label3;
        TEdit *Edit3;
        TButton *Button3;
        TGroupBox *GroupBox5;
        TLabel *Label4;
        TEdit *Edit4;
        TButton *Button11;
        TGroupBox *GroupBox4;
        TButton *Button6;
        TButton *Button7;
        TButton *Button8;
        TButton *Button9;
        TButton *Button10;
        TGroupBox *GroupBox2;
        TLabel *Label2;
        TEdit *Edit2;
        TButton *Button2;
        TButton *Button4;
        TButton *Button5;
        TBitBtn *BitBtn1;
        void __fastcall Button1Click(TObject *Sender);
        void __fastcall Button2Click(TObject *Sender);
        void __fastcall Button3Click(TObject *Sender);
        void __fastcall Button11Click(TObject *Sender);
        void __fastcall Button4Click(TObject *Sender);
        void __fastcall Button5Click(TObject *Sender);
        void __fastcall Button6Click(TObject *Sender);
        void __fastcall Button7Click(TObject *Sender);
        void __fastcall Button8Click(TObject *Sender);
        void __fastcall Button9Click(TObject *Sender);
        void __fastcall Button10Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
        __fastcall TForm1(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif

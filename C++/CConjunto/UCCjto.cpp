//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "UCCjto.h"
#include "CCjto.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
CConjunto<int,100> A,B;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Button1Click(TObject *Sender)
{
  if(Cjto->ItemIndex==0)
     A.Adicionar(StrToInt(Edit1->Text));
  else
     B.Adicionar(StrToInt(Edit1->Text));
  Edit1->Text = "";
  Edit1->SetFocus();
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button2Click(TObject *Sender)
{
  if(Cjto->ItemIndex==0)
    A.Eliminar(StrToInt(Edit2->Text));
  else
    B.Eliminar(StrToInt(Edit2->Text));

  Edit2->Text = "";
  Edit2->SetFocus();
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button3Click(TObject *Sender)
{
  if(Cjto->ItemIndex==0)
    ShowMessage(IntToStr(A.Obtener(StrToInt(Edit3->Text))));
  else
    ShowMessage(IntToStr(B.Obtener(StrToInt(Edit3->Text))));

  Edit3->Text = "";
  Edit3->SetFocus();
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button11Click(TObject *Sender)
{
   if(Cjto->ItemIndex==0){
      if(A.Pertenece(StrToInt(Edit4->Text)))
        ShowMessage("El elemento se encuentra en el Cjto");
      else
        ShowMessage("El elemento no se encuentra en el Cjto");
   }
   else{
      if(B.Pertenece(StrToInt(Edit4->Text)))
        ShowMessage("El elemento se encuentra en el Cjto");
      else
        ShowMessage("El elemento no se encuentra en el Cjto");
   }
   Edit4->Text = "";
   Edit4->SetFocus();
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button4Click(TObject *Sender)
{
    AnsiString CjtoA="{";
    for(int i=1 ;i<=A.Long();i++){
      CjtoA += IntToStr(A.Obtener(i));
      CjtoA += ",";
    }
    CjtoA+="}";
    ShowMessage(CjtoA);
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button5Click(TObject *Sender)
{
    AnsiString CjtoB="{";
    for(int i=1;i<=B.Long();i++){
      CjtoB += IntToStr(B.Obtener(i));
      CjtoB +=",";
    }
    CjtoB+="}";
    ShowMessage(CjtoB);
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button6Click(TObject *Sender)
{
  CConjunto<int, 100> C;

  C = A.Interseccion(B);

  AnsiString CjtoC="{";
  for(int i=1;i<=C.Long();i++){
    CjtoC += IntToStr(C.Obtener(i));
    CjtoC +=",";
  }
  CjtoC+="}";
  ShowMessage(CjtoC);
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button7Click(TObject *Sender)
{
  CConjunto<int, 100> C;

  C = A.Union(B);

  AnsiString CjtoC="{";
  for(int i=1;i<=C.Long();i++){
    CjtoC += IntToStr(C.Obtener(i));
    CjtoC +=",";
  }
  CjtoC+="}";
  ShowMessage(CjtoC);
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button8Click(TObject *Sender)
{
  if(A.Subconjunto(B)) ShowMessage("A es subcnjunto de B");
  else
   ShowMessage("A es subcnjunto de B");
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button9Click(TObject *Sender)
{
    CConjunto<int, 100> C;

  C = A.Diferencia(B);

  AnsiString CjtoC="{";
  for(int i=1;i<=C.Long();i++){
    CjtoC += IntToStr(C.Obtener(i));
    CjtoC +=",";
  }
  CjtoC+="}";
  ShowMessage(CjtoC);
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button10Click(TObject *Sender)
{
  if(A.Igual(B)) ShowMessage("A es Igual al B");
  else
   ShowMessage("A no es Igual al B");
}
//---------------------------------------------------------------------------

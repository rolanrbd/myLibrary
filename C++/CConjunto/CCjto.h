template<class T, int MaxLong=100>
class CConjunto{

 protected:
   T Items[MaxLong];
   int cant;

 public:
   CConjunto(){cant=0;};
   int Long(){return cant;};
   bool Pertenece(T x);
   void Adicionar(T x);
   void Eliminar(T x);
   T Obtener(int x);
   CConjunto Union(CConjunto B);
   CConjunto Diferencia(CConjunto B);
   CConjunto Interseccion(CConjunto B);
   bool Igual(CConjunto B);
   bool Subconjunto(CConjunto B);
   ~CConjunto(){cant = 0;};
};


template<class T, int MaxLong>
T CConjunto<T,MaxLong>::Obtener(int x){
 if(x>0 && x<cant)
   return Items[x-1];
}

template<class T, int MaxLong>
bool CConjunto<T,MaxLong>::Pertenece(T x){
 for(int i=0; i<cant; i++)
  if(Items[i]==x)
   return true;
 return false;
}

template<class T, int MaxLong>
void CConjunto<T,MaxLong>::Adicionar(T x){
 if(! Pertenece(x))
   Items[cant++]=x;
}

template<class T, int MaxLong>
void CConjunto<T,MaxLong>::Eliminar(T x){
int i=0;
 for(; i<cant; i++)
   if(Items[i]==x)break;

 for(; i<cant; i++)
   Items[i]=Items[i+1];
 cant--;
}

template<class T, int MaxLong>
CConjunto<T,MaxLong> CConjunto<T,MaxLong>::Union(CConjunto B){
CConjunto<T,MaxLong> C;

 for(int i; i<cant; i++)
  C.Adicionar(Items[i]);

  for(int i=1; i<=B.Long(); i++)
   C.Adicionar(B.Obtener(i));
  return C;
}

template<class T, int MaxLong>
CConjunto<T,MaxLong> CConjunto<T,MaxLong>::Diferencia(CConjunto B){

  CConjunto<T,MaxLong> D;

  for(int i=0; i<cant; i++)
   if(!B.Pertenece(Items[i]))
    D.Adicionar(Items[i]);

  return D;
}

template<class T, int MaxLong>
CConjunto<T,MaxLong> CConjunto<T,MaxLong>::Interseccion(CConjunto B){

  CConjunto<T,MaxLong> I;

  for(int i=0; i<cant; i++)
    if(B.Pertenece(Items[i]))
      I.Adicionar(Items[i]);

 return I;
}

template<class T, int MaxLong>
bool CConjunto<T,MaxLong>::Igual(CConjunto B){
  for(int i=0;i<cant;i++)
    if(B.Pertenece(Items[i]));
    else
      return false;

  return true;
}

template<class T, int MaxLong>
bool CConjunto<T,MaxLong>::Subconjunto(CConjunto B){
 int cont = 0;
 for(int i=0;i<cant;i++)
   if(B.Pertenece(Items[i]))
     cont++;

 if(cont <=cant-1)
   return true;
 else
   return false;
}

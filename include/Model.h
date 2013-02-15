
#ifndef MODEL_H
#define MODEL_H


class Model
{
    public:
    
        Model();
        ~Model();
        
        virtual void predict(RefArrayXd predictions, RefArrayXd parameters) = 0;
        
    protected:
    
    public:
    
};


#endif

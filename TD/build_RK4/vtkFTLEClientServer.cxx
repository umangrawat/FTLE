// ClientServer wrapper for vtkFTLE object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkFTLE.h"
#include "vtkSystemIncludes.h"
#include "vtkClientServerInterpreter.h"
#include "vtkClientServerStream.h"


vtkObjectBase *vtkFTLEClientServerNewCommand(void* /*ctx*/)
{
  return vtkFTLE::New();
}


int VTK_EXPORT vtkFTLECommand(vtkClientServerInterpreter *arlu, vtkObjectBase *ob, const char *method, const vtkClientServerStream& msg, vtkClientServerStream& resultStream, void* /*ctx*/)
{
  vtkFTLE *op = vtkFTLE::SafeDownCast(ob);
  if(!op)
    {
    vtkOStrStreamWrapper vtkmsg;
    vtkmsg << "Cannot cast " << ob->GetClassName() << " object to vtkFTLE.  "
           << "This probably means the class specifies the incorrect superclass in vtkTypeMacro.";
    resultStream.Reset();
    resultStream << vtkClientServerStream::Error
                 << vtkmsg.str() << 0 << vtkClientServerStream::End;
    return 0;
    }
  (void)arlu;
  if (!strcmp("New",method) && msg.GetNumberOfArguments(0) == 2)
    {
    vtkFTLE  *temp20;
      {
      temp20 = (op)->New();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("GetClassName",method) && msg.GetNumberOfArguments(0) == 2)
    {
    const char    *temp20;
      {
      temp20 = (op)->GetClassName();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("IsA",method) && msg.GetNumberOfArguments(0) == 3)
    {
    char    *temp0;
    int      temp20;
    if(msg.GetArgument(0, 2, &temp0))
      {
      temp20 = (op)->IsA(temp0);
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("NewInstance",method) && msg.GetNumberOfArguments(0) == 2)
    {
    vtkFTLE  *temp20;
      {
      temp20 = (op)->NewInstance();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SafeDownCast",method) && msg.GetNumberOfArguments(0) == 3)
    {
    vtkObject  *temp0;
    vtkFTLE  *temp20;
    if(vtkClientServerStreamGetArgumentObject(msg, 0, 2, &temp0, "vtkObject"))
      {
      temp20 = (op)->SafeDownCast(temp0);
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SettauInitial",method) && msg.GetNumberOfArguments(0) == 3)
    {
    double   temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SettauInitial(temp0);
      return 1;
      }
    }
  if (!strcmp("GettauInitial",method) && msg.GetNumberOfArguments(0) == 2)
    {
    double   temp20;
      {
      temp20 = (op)->GettauInitial();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SettauFinal",method) && msg.GetNumberOfArguments(0) == 3)
    {
    double   temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SettauFinal(temp0);
      return 1;
      }
    }
  if (!strcmp("GettauFinal",method) && msg.GetNumberOfArguments(0) == 2)
    {
    double   temp20;
      {
      temp20 = (op)->GettauFinal();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetmaxNumberOfSteps",method) && msg.GetNumberOfArguments(0) == 3)
    {
    int      temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SetmaxNumberOfSteps(temp0);
      return 1;
      }
    }
  if (!strcmp("GetmaxNumberOfSteps",method) && msg.GetNumberOfArguments(0) == 2)
    {
    int      temp20;
      {
      temp20 = (op)->GetmaxNumberOfSteps();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetintegrationStepSize",method) && msg.GetNumberOfArguments(0) == 3)
    {
    double   temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SetintegrationStepSize(temp0);
      return 1;
      }
    }
  if (!strcmp("GetintegrationStepSize",method) && msg.GetNumberOfArguments(0) == 2)
    {
    double   temp20;
      {
      temp20 = (op)->GetintegrationStepSize();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetoriginSeedGrid",method) && msg.GetNumberOfArguments(0) == 5)
    {
    double   temp0;
    double   temp1;
    double   temp2;
    if(msg.GetArgument(0, 2, &temp0) &&
      msg.GetArgument(0, 3, &temp1) &&
      msg.GetArgument(0, 4, &temp2))
      {
      op->SetoriginSeedGrid(temp0,temp1,temp2);
      return 1;
      }
    }
  if (!strcmp("SetoriginSeedGrid",method) && msg.GetNumberOfArguments(0) == 3)
    {
    double temp0[3];
    if(msg.GetArgument(0, 2, temp0, 3))
      {
      op->SetoriginSeedGrid(temp0);
      return 1;
      }
    }
  if (!strcmp("GetoriginSeedGrid",method) && msg.GetNumberOfArguments(0) == 2)
    {
    double  *temp20;
      {
      temp20 = (op)->GetoriginSeedGrid();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << vtkClientServerStream::InsertArray(temp20,3) << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetboundsSeedGrid",method) && msg.GetNumberOfArguments(0) == 5)
    {
    double   temp0;
    double   temp1;
    double   temp2;
    if(msg.GetArgument(0, 2, &temp0) &&
      msg.GetArgument(0, 3, &temp1) &&
      msg.GetArgument(0, 4, &temp2))
      {
      op->SetboundsSeedGrid(temp0,temp1,temp2);
      return 1;
      }
    }
  if (!strcmp("SetboundsSeedGrid",method) && msg.GetNumberOfArguments(0) == 3)
    {
    double temp0[3];
    if(msg.GetArgument(0, 2, temp0, 3))
      {
      op->SetboundsSeedGrid(temp0);
      return 1;
      }
    }
  if (!strcmp("GetboundsSeedGrid",method) && msg.GetNumberOfArguments(0) == 2)
    {
    double  *temp20;
      {
      temp20 = (op)->GetboundsSeedGrid();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << vtkClientServerStream::InsertArray(temp20,3) << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetdimensionSeedGrid",method) && msg.GetNumberOfArguments(0) == 5)
    {
    double   temp0;
    double   temp1;
    double   temp2;
    if(msg.GetArgument(0, 2, &temp0) &&
      msg.GetArgument(0, 3, &temp1) &&
      msg.GetArgument(0, 4, &temp2))
      {
      op->SetdimensionSeedGrid(temp0,temp1,temp2);
      return 1;
      }
    }
  if (!strcmp("SetdimensionSeedGrid",method) && msg.GetNumberOfArguments(0) == 3)
    {
    double temp0[3];
    if(msg.GetArgument(0, 2, temp0, 3))
      {
      op->SetdimensionSeedGrid(temp0);
      return 1;
      }
    }
  if (!strcmp("GetdimensionSeedGrid",method) && msg.GetNumberOfArguments(0) == 2)
    {
    double  *temp20;
      {
      temp20 = (op)->GetdimensionSeedGrid();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << vtkClientServerStream::InsertArray(temp20,3) << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetuseGPU",method) && msg.GetNumberOfArguments(0) == 3)
    {
    bool   temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SetuseGPU(temp0);
      return 1;
      }
    }
  if (!strcmp("GetuseGPU",method) && msg.GetNumberOfArguments(0) == 2)
    {
    bool   temp20;
      {
      temp20 = (op)->GetuseGPU();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetshowStreamLines",method) && msg.GetNumberOfArguments(0) == 3)
    {
    bool   temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SetshowStreamLines(temp0);
      return 1;
      }
    }
  if (!strcmp("GetshowStreamLines",method) && msg.GetNumberOfArguments(0) == 2)
    {
    bool   temp20;
      {
      temp20 = (op)->GetshowStreamLines();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("printSeedGridProperties",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->printSeedGridProperties();
      return 1;
      }
    }

  {
    const char* commandName = "vtkMultiTimeStepAlgorithm";
    if (arlu->HasCommandFunction(commandName) &&
        arlu->CallCommandFunction(commandName, op, method, msg, resultStream)) { return 1; }
  }
  if(resultStream.GetNumberOfMessages() > 0 &&
     resultStream.GetCommand(0) == vtkClientServerStream::Error &&
     resultStream.GetNumberOfArguments(0) > 1)
    {
    /* A superclass wrapper prepared a special message. */
    return 0;
    }
  vtkOStrStreamWrapper vtkmsg;
  vtkmsg << "Object type: vtkFTLE, could not find requested method: \""
         << method << "\"\nor the method was called with incorrect arguments.\n";
  resultStream.Reset();
  resultStream << vtkClientServerStream::Error
               << vtkmsg.str() << vtkClientServerStream::End;
  vtkmsg.rdbuf()->freeze(0);
  return 0;
}


//-------------------------------------------------------------------------auto
void VTK_EXPORT vtkFTLE_Init(vtkClientServerInterpreter* csi)
{
  static vtkClientServerInterpreter* last = NULL;
  if(last != csi)
    {
    last = csi;
    csi->AddNewInstanceFunction("vtkFTLE", vtkFTLEClientServerNewCommand);
    csi->AddCommandFunction("vtkFTLE", vtkFTLECommand);
    }
}

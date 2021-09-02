#include "suids.hpp"

namespace emodlib
{
    namespace suids
    {
        distributed_generator::distributed_generator( int _rank, int _numtasks )
        : rank( _rank )
        , numtasks( _numtasks )
        {
            next_suid.data = rank + 1; // +1 ensures that NIL will never be generated
        }

        suid distributed_generator::operator()()
        {
            suid tmp = next_suid;
            next_suid.data += numtasks;
            return tmp;
        }

    }
}

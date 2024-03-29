#ifndef OCCA_INTERNAL_LANG_TYPE_UNION_HEADER
#define OCCA_INTERNAL_LANG_TYPE_UNION_HEADER

#include <occa/internal/lang/type/type.hpp>

namespace occa {
  namespace lang {
    class union_t : public type_t {
    public:
      variableVector fields;

      union_t();
      union_t(identifierToken &nameToken);

      union_t(const union_t &other);

      virtual int type() const;
      virtual type_t& clone() const;

      virtual dtype_t dtype() const;

      void addField(variable_t &var);
      void addFields(variableVector &fields_);

      virtual void printDeclaration(printer &pout) const;
    };
  }
}

#endif

#!/bin/sh
#
# Edit encoding-sided-in.c and test-sided-in.c
# Run this thing.
# Compile in a Cygwin32 window with:
#	clang test.c 2>&1 | less  -R
#

cp encoding-sided-in.c encoding-sided.c
cp test-sided-in.c test-sided.c
cp test-mul-in.c test-mul.c
fn="encoding-sided.c test-sided.c test-mul.c"
sed -i 's/encoding->new_elt(/encoding_new_elt(/g' $fn
sed -i 's/encoding->del_elt(/encoding_del_elt(/g' $fn
sed -i 's/encoding->encode(/encoding_encode(/g' $fn
sed -i 's/encoding->equal(/encoding_equal(/g' $fn
sed -i 's/encoding->add(/encoding_add(/g' $fn
sed -i 's/encoding->sub(/encoding_sub(/g' $fn
sed -i 's/encoding->new_prod(/encoding_new_prod(encoding/g' $fn
sed -i 's/encoding->mul(/encoding_mul(/g' $fn
sed -i 's/field->newElt(/field_newElt(field/g' $fn
sed -i 's/field->newEltArray(/field_newEltArray(/g' $fn
sed -i 's/field->assignRandomElt(/field_assignRandomElt(/g' $fn
sed -i 's/field->delElt(/field_delElt(/g' $fn
sed -i 's/field->one(/field_one(/g' $fn
sed -i 's/field->add(/field_add(/g' $fn
sed -i 's/field->sub(/field_sub(/g' $fn
sed -i 's/field->mul(/field_mul(/g' $fn

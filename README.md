# scdata

A re-implementation of my Rustody::singlecelldata::SingleCellData class.
In comparison to that this class can now collect per gene+UMI f32 values and stores them as mean f32 in the saved files.
In order to get this you need to create the object using the ``MatrixValueType::Real`` like that.

```
use scdata::Scdata;
use scdata::cell_data::GeneUmiHash;
use crate::cell_data::GeneUmiHash;

let mut celldata = Scdata::new( 1, MatrixValueType::Real ); // only one thread here

for i in 0..20 {
    celldata.try_insert( &13452355_u64, GeneUmiHash(2, umi as u64), i as f32, &mut report)
}
```

This will get the final value of 9.5 in the saved matrix marked table.

If you want to also save the total read count of the data you need to modify the celldata.value_type( MatrixValueType::Integer) and re-export the data.



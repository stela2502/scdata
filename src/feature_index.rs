// feature_index.rs

pub trait FeatureIndex: Send + Sync {
    /// name to id translation
    fn feature_name(&self, feature_id: u64) -> &str;
    /// id to name translation
    fn feature_id(&self, name: &str) -> Option<u64>;
    /// one line in the features.tsv file
    fn to_10x_feature_line(&self, feature_id: u64) -> String;
    /// Feature ids in canonical export order.
    fn ordered_feature_ids(&self) -> Vec<u64>;
}
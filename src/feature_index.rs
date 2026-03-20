// feature_index.rs

pub trait FeatureIndex: Send + Sync {
    fn feature_name(&self, feature_id: u64) -> &str;
    fn to_10x_feature_line(&self, feature_id: u64) -> String;
}
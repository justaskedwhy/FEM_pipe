// cpp/include/element_registry.hpp
#pragma once
#include "element_kernel.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <functional>
#include <vector>

namespace fem {

class ElementRegistry {
public:
    using Factory = std::function<std::unique_ptr<ElementKernel>()>;

    static ElementRegistry& instance() {
        static ElementRegistry instance_;
        return instance_;
    }

    void register_kernel(const std::string& abaqus_name, Factory f) {
        factories_[abaqus_name] = f;
    }

    std::unique_ptr<ElementKernel> create(const std::string& abaqus_name) const {
        auto it = factories_.find(abaqus_name);
        if (it == factories_.end()) {
            throw std::runtime_error("[ElementRegistry] Unregistered element: " + abaqus_name);
        }
        return it->second();
    }

    bool has(const std::string& abaqus_name) const {
        return factories_.find(abaqus_name) != factories_.end();
    }

    std::vector<std::string> registered_names() const {
        std::vector<std::string> names;
        for (const auto& [name, _] : factories_) {
            names.push_back(name);
        }
        return names;
    }

private:
    std::unordered_map<std::string, Factory> factories_;
};

#define CAT_(a, b) a ## b
#define CAT(a, b) CAT_(a, b)

#define REGISTER_ELEMENT(ABAQUS_NAME, CLASS_NAME)                              \
    REGISTER_ELEMENT_IMPL(ABAQUS_NAME, CLASS_NAME, CAT(_registered_, __LINE__))

#define REGISTER_ELEMENT_IMPL(ABAQUS_NAME, CLASS_NAME, VAR)                    \
    namespace {                                                                \
        const bool VAR = []() {                                                \
            fem::ElementRegistry::instance().register_kernel(                  \
                ABAQUS_NAME,                                                   \
                []() {                                                         \
                    auto k = std::make_unique<CLASS_NAME>();                   \
                    k->set_registry_name(ABAQUS_NAME);                         \
                    return k;                                                  \
                });                                                            \
            return true;                                                       \
        }();                                                                   \
    }

} // namespace fem